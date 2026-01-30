"""
RAVUB Calculator - Ranked AV Upper Bound Analysis Tool

This script processes ranked ballot data and computes RAVUB (Ranked AV Upper Bound)
for each major candidate in an election, outputting results as JSON and generating
dark-mode visualizations.
"""

import json
import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import Optional
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


class BallotProcessor:
    """Processes ranked ballots and handles edge cases like overvotes, undervotes, and duplicates."""
    
    def __init__(self, major_candidates: list[str], write_in_markers: list[str]):
        self.major_candidates = set(major_candidates)
        self.major_candidates_list = major_candidates  # Preserve order
        self.write_in_markers = set(write_in_markers)
        self.stats = {
            "total_ballots": 0,
            "valid_ballots": 0,
            "empty_ballots": 0,
            "overvotes_corrected": 0,
            "undervotes": 0,
            "write_ins_encountered": 0
        }
    
    def is_write_in(self, candidate: str) -> bool:
        """Check if a candidate is a write-in."""
        return candidate in self.write_in_markers or candidate.lower().startswith("write")
    
    def clean_ballot(self, rankings: list[str]) -> list[str]:
        """
        Clean a ballot by:
        1. Removing empty entries
        2. Removing duplicate candidates (overvote correction)
        3. Keeping track of statistics
        
        Returns the cleaned list of rankings.
        """
        cleaned = []
        seen = set()
        has_overvote = False
        
        for rank in rankings:
            # Skip empty entries
            if not rank or rank.strip() == "":
                continue
            
            candidate = rank.strip()
            
            # Skip if we've already seen this candidate (overvote)
            if candidate in seen:
                has_overvote = True
                continue
            
            seen.add(candidate)
            cleaned.append(candidate)
            
            # Track write-ins
            if self.is_write_in(candidate):
                self.stats["write_ins_encountered"] += 1
        
        if has_overvote:
            self.stats["overvotes_corrected"] += 1
        
        return cleaned
    
    def extract_major_candidate_order(self, rankings: list[str]) -> list[str]:
        """
        Extract the implied ordering of major candidates from a ballot.
        Non-major candidates (including write-ins) are filtered out.
        
        Returns list of major candidates in the order they appear on the ballot.
        """
        major_order = []
        for candidate in rankings:
            if candidate in self.major_candidates:
                major_order.append(candidate)
        return major_order
    
    def get_implicit_ranking(self, rankings: list[str]) -> dict[str, int]:
        """
        Get the implicit ranking of all major candidates.
        Candidates not on the ballot are ranked last (after all ranked majors).
        
        Returns dict mapping candidate -> rank (1-indexed, lower is better)
        """
        major_order = self.extract_major_candidate_order(rankings)
        ranking = {}
        
        # Assign ranks to candidates that appear on the ballot
        for i, candidate in enumerate(major_order):
            ranking[candidate] = i + 1
        
        # Assign last rank to unranked major candidates
        last_rank = len(major_order) + 1
        for candidate in self.major_candidates:
            if candidate not in ranking:
                ranking[candidate] = last_rank
        
        return ranking
    
    def is_implicitly_ranked_above(self, rankings: list[str], candidate_x: str, candidate_y: str) -> bool:
        """
        Check if candidate X is implicitly ranked above candidate Y.
        
        X is implicitly ranked above Y if:
        1. Both are ranked and X appears before Y, OR
        2. X is ranked and Y is not ranked
        """
        implicit = self.get_implicit_ranking(rankings)
        return implicit.get(candidate_x, float('inf')) < implicit.get(candidate_y, float('inf'))
    
    def candidate_is_ranked(self, rankings: list[str], candidate: str) -> bool:
        """Check if a candidate appears on the ballot."""
        return candidate in rankings


class RAVUBCalculator:
    """Calculates Ranked AV Upper Bound for elections."""
    
    def __init__(self, config_path: str = "elections_config.json"):
        self.config_path = config_path
        self.config = self._load_config()
        self.results = {}
    
    def _load_config(self) -> dict:
        """Load election configuration from JSON file."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def _load_ballots(self, election_id: str) -> tuple[list[tuple[list[str], int]], BallotProcessor]:
        """
        Load and process ballots from CSV file.
        
        Returns tuple of (processed_ballots, processor) where processed_ballots
        is a list of (cleaned_rankings, count) tuples.
        """
        election_config = self.config["elections"][election_id]
        data_file = Path("data") / election_config["data_file"]
        
        processor = BallotProcessor(
            major_candidates=election_config["major_candidates"],
            write_in_markers=election_config.get("write_in_markers", ["Write-in"])
        )
        
        ballots = []
        
        with open(data_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                # Extract count
                count = int(row.get("count", 1))
                processor.stats["total_ballots"] += count
                
                # Extract rankings (handle variable number of rank columns)
                rankings = []
                for key in sorted(row.keys()):
                    if key.startswith("rank") and key != "count":
                        value = row[key]
                        if value:
                            rankings.append(value)
                
                # Clean the ballot
                cleaned = processor.clean_ballot(rankings)
                
                if not cleaned:
                    processor.stats["empty_ballots"] += count
                else:
                    processor.stats["valid_ballots"] += count
                    ballots.append((cleaned, count))
        
        return ballots, processor
    
    def compute_ravub(self, election_id: str, profile_candidate: str, 
                      ballots: list[tuple[list[str], int]], 
                      processor: BallotProcessor) -> dict:
        """
        Compute RAVUB (Ranked AV Upper Bound) for a given candidate.
        
        RAVUB rule:
        1. If candidate X is ranked on a voter's ballot, and X is implicitly ranked 
           above some major candidate Y in M-{X}, then approve X and all candidates 
           ranked above X.
        2. Otherwise, only approve the top-ranked candidate (even if it's a write-in,
           meaning no major candidate gets approved from this ballot).
        """
        major_candidates = processor.major_candidates_list
        approval_counts = defaultdict(int)
        
        for rankings, count in ballots:
            if not rankings:
                continue
            
            top_candidate = rankings[0]
            approvals = set()
            
            # Check if profile candidate is ranked and implicitly above some other major
            profile_ranked = processor.candidate_is_ranked(rankings, profile_candidate)
            
            if profile_ranked:
                # Check if profile candidate is implicitly ranked above any other major
                implicitly_above_some_major = False
                for other_major in major_candidates:
                    if other_major != profile_candidate:
                        if processor.is_implicitly_ranked_above(rankings, profile_candidate, other_major):
                            implicitly_above_some_major = True
                            break
                
                if implicitly_above_some_major:
                    # Approve profile candidate and all major candidates ranked above them
                    major_order = processor.extract_major_candidate_order(rankings)
                    
                    if profile_candidate in major_order:
                        profile_pos = major_order.index(profile_candidate)
                        # Approve profile candidate and all majors above
                        for i in range(profile_pos + 1):
                            approvals.add(major_order[i])
                else:
                    # Profile candidate is ranked last among majors, only approve top candidate
                    if top_candidate in processor.major_candidates:
                        approvals.add(top_candidate)
            else:
                # Profile candidate not ranked, only approve top candidate
                if top_candidate in processor.major_candidates:
                    approvals.add(top_candidate)
            
            # Tally approvals (only count major candidates)
            for candidate in approvals:
                if candidate in processor.major_candidates:
                    approval_counts[candidate] += count
        
        return dict(approval_counts)
    
    def compute_first_choice_counts(self, ballots: list[tuple[list[str], int]], 
                                       processor: BallotProcessor) -> dict:
        """
        Compute first-choice (plurality) vote counts for each major candidate.
        
        Only counts ballots where a major candidate is literally ranked first.
        Ballots with write-ins or non-major candidates ranked first do not
        contribute to any major candidate's first-choice count.
        
        Returns dict mapping candidate -> first-choice vote count
        """
        first_choice_counts = defaultdict(int)
        
        for rankings, count in ballots:
            if not rankings:
                continue
            
            # Only count if rank 1 is a major candidate
            top_candidate = rankings[0]
            if top_candidate in processor.major_candidates:
                first_choice_counts[top_candidate] += count
        
        return dict(first_choice_counts)
    
    def analyze_election(self, election_id: str) -> dict:
        """
        Analyze an election and compute RAVUBs for all major candidates.
        
        Returns complete analysis with all statistics.
        """
        election_config = self.config["elections"][election_id]
        ballots, processor = self._load_ballots(election_id)
        
        major_candidates = election_config["major_candidates"]
        total_ballots = processor.stats["total_ballots"]
        
        # Compute first-choice counts
        first_choice_counts = self.compute_first_choice_counts(ballots, processor)
        
        results = {
            "election_id": election_id,
            "election_name": election_config["name"],
            "date": election_config["date"],
            "description": election_config.get("description", ""),
            "method": "Ranked AV Upper Bound (RAVUB)",
            "total_ballots": total_ballots,
            "valid_ballots": processor.stats["valid_ballots"],
            "empty_ballots": processor.stats["empty_ballots"],
            "irv_winner": election_config.get("irv_winner", "Unknown"),
            "ballot_stats": processor.stats,
            "major_candidates": major_candidates,
            "first_choice_counts": {
                candidate: {
                    "votes": first_choice_counts.get(candidate, 0),
                    "percentage": round(first_choice_counts.get(candidate, 0) / total_ballots * 100, 2) if total_ballots > 0 else 0
                }
                for candidate in major_candidates
            },
            "profiles": {}
        }
        
        for profile_candidate in major_candidates:
            approval_counts = self.compute_ravub(election_id, profile_candidate, ballots, processor)
            
            # Sort by votes descending
            sorted_results = sorted(
                [(c, approval_counts.get(c, 0)) for c in major_candidates],
                key=lambda x: (-x[1], x[0])
            )
            
            # Determine winner and runner-up
            winner = sorted_results[0][0]
            winner_votes = sorted_results[0][1]
            runner_up = sorted_results[1][0] if len(sorted_results) > 1 else None
            runner_up_votes = sorted_results[1][1] if len(sorted_results) > 1 else 0
            
            wins_own_profile = (winner == profile_candidate)
            profile_votes = approval_counts.get(profile_candidate, 0)
            
            # Calculate margins
            if wins_own_profile:
                margin = winner_votes - runner_up_votes
                margin_opponent = runner_up
            else:
                margin = profile_votes - winner_votes  # Will be negative
                margin_opponent = winner
            
            margin_percentage = (margin / total_ballots * 100) if total_ballots > 0 else 0
            
            profile_result = {
                "profile_candidate": profile_candidate,
                "winner": winner,
                "wins_own_profile": wins_own_profile,
                "margin": margin,
                "margin_percentage": round(margin_percentage, 2),
                "margin_opponent": margin_opponent,
                "approval_votes": [
                    {
                        "candidate": candidate,
                        "votes": votes,
                        "percentage": round(votes / total_ballots * 100, 2) if total_ballots > 0 else 0
                    }
                    for candidate, votes in sorted_results
                ]
            }
            
            results["profiles"][profile_candidate] = profile_result
        
        return results
    
    def analyze_all_elections(self) -> dict:
        """Analyze all elections defined in the config."""
        all_results = {}
        
        for election_id in self.config["elections"]:
            print(f"Analyzing {election_id}...")
            all_results[election_id] = self.analyze_election(election_id)
        
        return all_results
    
    def save_results(self, results: dict, output_dir: str = "output/UBs"):
        """Save results to JSON files."""
        os.makedirs(output_dir, exist_ok=True)
        
        for election_id, election_results in results.items():
            output_path = Path(output_dir) / f"{election_id}_ravub.json"
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(election_results, f, indent=2, ensure_ascii=False)
            print(f"Saved: {output_path}")
        
        # Also save combined results
        combined_path = Path(output_dir) / "all_ravub_results.json"
        with open(combined_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"Saved: {combined_path}")


class RAVUBVisualizer:
    """Generates dark-mode visualizations for RAVUB results."""
    
    # Dark mode color scheme
    COLORS = {
        'background': '#1a1a2e',
        'surface': '#16213e',
        'text': '#eaeaea',
        'text_secondary': '#a0a0a0',
        'grid': '#2a2a4a',
        'winner': '#4ade80',
        'loser': '#f87171',
        'neutral': '#60a5fa',
        'bar_colors': ['#818cf8', '#34d399', '#fbbf24', '#f472b6', '#22d3d8', '#a78bfa', '#fb923c']
    }
    
    def __init__(self, output_dir: str = "output/UBs"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set dark mode style
        plt.style.use('dark_background')
        plt.rcParams.update({
            'figure.facecolor': self.COLORS['background'],
            'axes.facecolor': self.COLORS['surface'],
            'axes.edgecolor': self.COLORS['grid'],
            'axes.labelcolor': self.COLORS['text'],
            'text.color': self.COLORS['text'],
            'xtick.color': self.COLORS['text'],
            'ytick.color': self.COLORS['text'],
            'grid.color': self.COLORS['grid'],
            'font.family': 'sans-serif',
            'font.size': 11
        })
    
    def generate_profile_chart(self, election_results: dict, profile_candidate: str) -> str:
        """Generate a bar chart for a single RAVUB profile."""
        profile = election_results["profiles"][profile_candidate]
        election_name = election_results["election_name"]
        total_ballots = election_results["total_ballots"]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor(self.COLORS['background'])
        ax.set_facecolor(self.COLORS['surface'])
        
        # Extract data
        candidates = [v["candidate"] for v in profile["approval_votes"]]
        votes = [v["votes"] for v in profile["approval_votes"]]
        percentages = [v["percentage"] for v in profile["approval_votes"]]
        
        # Create bar colors
        bar_colors = []
        for i, candidate in enumerate(candidates):
            if candidate == profile["winner"]:
                bar_colors.append(self.COLORS['winner'])
            elif candidate == profile_candidate and not profile["wins_own_profile"]:
                bar_colors.append(self.COLORS['loser'])
            else:
                bar_colors.append(self.COLORS['bar_colors'][i % len(self.COLORS['bar_colors'])])
        
        # Create horizontal bar chart
        y_pos = np.arange(len(candidates))
        bars = ax.barh(y_pos, votes, color=bar_colors, edgecolor='white', linewidth=0.5)
        
        # Add vote counts and percentages on bars
        for i, (bar, vote, pct) in enumerate(zip(bars, votes, percentages)):
            width = bar.get_width()
            label = f'{vote:,} ({pct:.1f}%)'
            
            # Position text inside or outside bar depending on size
            if width > max(votes) * 0.3:
                ax.text(width - max(votes) * 0.02, bar.get_y() + bar.get_height()/2,
                       label, ha='right', va='center', color='white', fontweight='bold', fontsize=10)
            else:
                ax.text(width + max(votes) * 0.01, bar.get_y() + bar.get_height()/2,
                       label, ha='left', va='center', color=self.COLORS['text'], fontsize=10)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(candidates)
        ax.invert_yaxis()
        
        ax.set_xlabel('Approval Votes', fontsize=12)
        ax.set_title(f'RAVUB Profile: {profile_candidate}', fontsize=14, fontweight='bold', pad=15)
        
        # Add subtitle with election info
        subtitle = f'{election_name} | Total Ballots: {total_ballots:,}'
        ax.text(0.5, 1.02, subtitle, transform=ax.transAxes, ha='center', 
                fontsize=10, color=self.COLORS['text_secondary'])
        
        # Add result annotation
        if profile["wins_own_profile"]:
            result_text = f"✓ {profile_candidate} WINS their profile"
            result_color = self.COLORS['winner']
        else:
            result_text = f"✗ {profile_candidate} LOSES to {profile['winner']}"
            result_color = self.COLORS['loser']
        
        margin_text = f"Margin: {profile['margin']:+,} ({profile['margin_percentage']:+.2f}%)"
        
        # Add result box
        props = dict(boxstyle='round,pad=0.5', facecolor=self.COLORS['surface'], 
                     edgecolor=result_color, linewidth=2)
        ax.text(0.98, 0.02, f"{result_text}\n{margin_text}", transform=ax.transAxes,
               fontsize=10, va='bottom', ha='right', bbox=props, color=result_color)
        
        ax.set_xlim(0, max(votes) * 1.15)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # Save figure
        safe_candidate = profile_candidate.replace(" ", "_").replace(",", "").replace(".", "")
        filename = f"{election_results['election_id']}_{safe_candidate}_ravub.png"
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=150, bbox_inches='tight', facecolor=self.COLORS['background'])
        plt.close()
        
        return str(filepath)
    
    def generate_all_charts(self, all_results: dict) -> dict:
        """Generate charts for all elections and profiles."""
        chart_paths = {}
        
        for election_id, election_results in all_results.items():
            chart_paths[election_id] = {}
            
            for profile_candidate in election_results["profiles"]:
                filepath = self.generate_profile_chart(election_results, profile_candidate)
                chart_paths[election_id][profile_candidate] = filepath
                print(f"Generated: {filepath}")
        
        return chart_paths
    
    def generate_summary_chart(self, election_results: dict) -> str:
        """Generate a summary chart showing all profiles for an election."""
        election_id = election_results["election_id"]
        election_name = election_results["election_name"]
        profiles = election_results["profiles"]
        major_candidates = election_results["major_candidates"]
        
        fig, ax = plt.subplots(figsize=(12, 7))
        fig.patch.set_facecolor(self.COLORS['background'])
        ax.set_facecolor(self.COLORS['surface'])
        
        n_candidates = len(major_candidates)
        x = np.arange(n_candidates)
        width = 0.8 / n_candidates
        
        # For each profile candidate, show the votes each candidate gets
        for i, profile_candidate in enumerate(major_candidates):
            profile = profiles[profile_candidate]
            votes_dict = {v["candidate"]: v["percentage"] for v in profile["approval_votes"]}
            votes = [votes_dict.get(c, 0) for c in major_candidates]
            
            offset = (i - n_candidates/2 + 0.5) * width
            color = self.COLORS['bar_colors'][i % len(self.COLORS['bar_colors'])]
            bars = ax.bar(x + offset, votes, width, label=f'{profile_candidate} profile', 
                         color=color, edgecolor='white', linewidth=0.3)
        
        ax.set_xlabel('Candidates', fontsize=12)
        ax.set_ylabel('Approval Percentage (%)', fontsize=12)
        ax.set_title(f'RAVUB Summary: {election_name}', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(major_candidates, rotation=45, ha='right')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        filename = f"{election_id}_ravub_summary.png"
        filepath = self.output_dir / filename
        plt.savefig(filepath, dpi=150, bbox_inches='tight', facecolor=self.COLORS['background'])
        plt.close()
        
        return str(filepath)


def main():
    """Main entry point for RAVUB analysis."""
    print("=" * 60)
    print("RAVUB Calculator - Ranked AV Upper Bound Analysis")
    print("=" * 60)
    
    # Initialize calculator
    calculator = RAVUBCalculator()
    
    # Analyze all elections
    results = calculator.analyze_all_elections()
    
    # Save JSON results
    calculator.save_results(results)
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    visualizer = RAVUBVisualizer()
    
    for election_id, election_results in results.items():
        # Generate individual profile charts
        for profile_candidate in election_results["profiles"]:
            visualizer.generate_profile_chart(election_results, profile_candidate)
        
        # Generate summary chart
        visualizer.generate_summary_chart(election_results)
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)
    
    # Print summary
    for election_id, election_results in results.items():
        print(f"\n{election_results['election_name']}:")
        print(f"  Total ballots: {election_results['total_ballots']:,}")
        print(f"  IRV Winner: {election_results['irv_winner']}")
        print(f"  RAVUB Results:")
        for profile_candidate, profile in election_results["profiles"].items():
            status = "✓ WINS" if profile["wins_own_profile"] else f"✗ loses to {profile['winner']}"
            print(f"    {profile_candidate}: {status} (margin: {profile['margin']:+,})")


if __name__ == "__main__":
    main()
