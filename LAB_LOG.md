# Laboratory log

*The Effect of AI Companion Voice on Loneliness* — Registered Report, *Nature Human Behaviour*.

Stage 1 protocol (accepted in principle) registered at
<https://doi.org/10.6084/m9.figshare.29504894>, posted 2025-09-04. **All data collection for
both studies took place after in-principle acceptance and after that registration date**, as
the ranges below show.

Dates are taken from the `StartDate` field of the survey exports in this repository and can
be re-derived from them directly. Ranges are inclusive and given in the survey platform's
recorded local dates.

## Study 1

| | |
|---|---|
| Data collection | 2025-12-17 to 2026-01-20 |
| Platform | CloudResearch Connect |
| Session records collected | 4,743 |
| Final analysed sample | 1,260 (after attention, loneliness-screening and comprehension exclusions) |
| File | `e1/data.csv` |

Study 1 was a single-session experiment; each participant contributed one session.

## Study 2

Ten consecutive daily sessions per participant. Because participants started on staggered
dates and each day's survey had a completion window, the per-day ranges overlap.

| Day | File | Collection dates | Session records |
|---|---|---|---|
| 1 | `e2/d1.csv` | 2026-05-14 to 2026-05-18 | 3,147 |
| 2 | `e2/d2.csv` | 2026-05-14 to 2026-05-25 | 813 |
| 3 | `e2/d3.csv` | 2026-05-15 to 2026-05-25 | 675 |
| 4 | `e2/d4.csv` | 2026-05-16 to 2026-05-23 | 613 |
| 5 | `e2/d5.csv` | 2026-05-17 to 2026-05-25 | 579 |
| 6 | `e2/d6.csv` | 2026-05-18 to 2026-05-24 | 557 |
| 7 | `e2/d7.csv` | 2026-05-19 to 2026-05-24 | 539 |
| 8 | `e2/d8.csv` | 2026-05-20 to 2026-05-27 | 529 |
| 9 | `e2/d9.csv` | 2026-05-21 to 2026-05-27 | 492 |
| 10 | `e2/d10.csv` | 2026-05-22 to 2026-05-27 | 491 |

Overall Study 2 collection ran 2026-05-14 to 2026-05-27. The Day 1 record count includes
participants screened out at baseline for insufficient loneliness, who did not continue. The
final analysed sample is the 461 participants who completed all ten days (text = 153,
voice = 137, control = 171).

## Notes

Analyses were conducted after data collection closed for each study. The analysis scripts and
the exclusion chain that produces these final samples are in `e1/e1.R` and `e2/e2.R`; see
[README.md](README.md) for how to run them and for the conversation data that is withheld
under our IRB-approved protocol.
