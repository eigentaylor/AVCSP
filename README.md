# RAVUB Analysis

**Ranked AV Upper Bound (RAVUB)** — Deducing possible Approval Voting outcomes from ranked ballot data.

## Overview

This project analyzes real-world ranked-choice voting (RCV) elections to determine whether candidates could win under Approval Voting, using a theoretical construction called the **Ranked AV Upper Bound (RAVUB)**.

Based on Steven Brams' *AV Critical Strategy Profile* from [Mathematics and Democracy (2008)](https://press.princeton.edu/books/paperback/9780691133218/mathematics-and-democracy), the RAVUB provides an upper bound on each candidate's performance under Approval Voting derived from ranked ballot data.

## Key Concepts

### Implicit Ranking

A candidate X is **implicitly ranked** above candidate Y on a voter's ballot if:
1. Both X and Y are ranked, and X appears before Y, **OR**
2. X is ranked and Y is not ranked at all

### RAVUB Definition

Let **M** be the set of major candidates. The **Ranked AV Upper Bound** for candidate X is computed as follows:

1. If a voter ranks X and X is implicitly ranked above some other major candidate Y ∈ M\\{X}, then the voter approves X and all candidates ranked above X.
2. Otherwise, the voter only approves their top-ranked candidate.

This construction maximizes X's approval count under sincere voting strategies.

### Why "Major Candidates"?

We specify "major candidates" to treat ballots with the same implied preference structure identically. For example, analyzing the Alaska 2022 special election with major candidates {Peltola, Begich, Palin}, the following ballots all imply **Palin > Begich > Peltola**:

- Palin > Begich
- Palin > Begich > Peltola
- Palin > Begich > Write-in
- Palin > Write-in > Begich
- Palin > Begich > Write-in > Peltola

### What RAVUB Can Prove

- **If X does not win under their own RAVUB profile**, then X cannot win under any sincere Approval Voting scenario.
- **If X is a Condorcet winner**, X will also win under their RAVUB profile.
- **Even a Condorcet loser** may win under their own RAVUB profile (though they'd lose under others).

## Elections Analyzed

| Election | Date | IRV Winner |
|----------|------|------------|
| Alaska 2022 Special General | 2022-08-16 | Mary Peltola |
| Minneapolis Ward 2 City Council | 2021-11-02 | Robin Wonsley Worlobah |
| NYC 2025 Democratic Mayoral Primary | 2025-06-24 | Zohran Kwame Mamdani |

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/AVCSP.git
cd AVCSP

# Create virtual environment (optional but recommended)
python -m venv .venv
.venv\Scripts\activate  # Windows
# source .venv/bin/activate  # Linux/Mac

# Install dependencies
pip install matplotlib numpy
```

## Usage

### Running the Analysis

```bash
python ravub_calculator.py
```

This will:
1. Load ballot data from `data/` folder
2. Process each election using the config in `elections_config.json`
3. Compute RAVUB for all major candidates
4. Output JSON results to `output/UBs/`
5. Generate dark-mode visualization charts

### Viewing Results

Open `index.html` in a browser to view an interactive dashboard with:
- Election and candidate profile selection
- Approval vote counts and percentages
- Win/loss status with margins
- Visualization charts

For GitHub Pages: push to your repository and enable Pages in settings.

## Project Structure

```
AVCSP/
├── data/                          # Ballot data CSV files
│   ├── alaska_2022_ballots_named.csv
│   ├── mn_ward2_2021_ballots_named.csv
│   └── nyc_2025_mayor_ballots_named.csv
├── output/
│   └── UBs/                       # Generated results
│       ├── all_ravub_results.json
│       ├── *_ravub.json           # Per-election results
│       └── *.png                  # Visualization charts
├── elections_config.json          # Election configuration
├── ravub_calculator.py            # Main analysis script
├── index.html                     # Interactive results viewer
└── README.md
```

## Configuration

Edit `elections_config.json` to add new elections:

```json
{
  "elections": {
    "election_id": {
      "name": "Election Display Name",
      "date": "YYYY-MM-DD",
      "data_file": "filename.csv",
      "description": "Brief description",
      "major_candidates": ["Candidate A", "Candidate B"],
      "irv_winner": "Candidate A",
      "write_in_markers": ["Write-in"],
      "max_rankings": 5
    }
  }
}
```

### Ballot Data Format

CSV files should have columns: `rank1`, `rank2`, ..., `rankN`, `count`

```csv
rank1,rank2,rank3,count
"Candidate A","Candidate B",,1500
"Candidate B","Candidate A","Candidate C",800
```

## Handling Edge Cases

The calculator properly handles:
- **Overvotes**: Duplicate rankings (e.g., "Begich > Begich") are deduplicated
- **Undervotes**: Empty ballots are tracked but excluded from calculations
- **Write-ins**: Configurable markers; excluded from major candidate analysis
- **Partial ballots**: Unranked major candidates treated as implicitly ranked last

## Output Format

```json
{
  "election_id": "alaska_2022",
  "election_name": "Alaska 2022 Special General Election",
  "total_ballots": 192289,
  "profiles": {
    "Candidate Name": {
      "winner": "Winning Candidate",
      "wins_own_profile": true,
      "margin": 12345,
      "margin_percentage": 6.42,
      "approval_votes": [
        {"candidate": "X", "votes": 100000, "percentage": 52.0}
      ]
    }
  }
}
```

## References

- Brams, S. J. (2008). *Mathematics and Democracy*. Princeton University Press.
- [ranked.vote](https://ranked.vote) — RCV election data and analysis

## License

MIT License