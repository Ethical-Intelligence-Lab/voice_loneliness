# The Effect of AI Companion Voice on Loneliness

Analysis code and de-identified data for the Registered Report *The Effect of AI Companion Voice on Loneliness* (Stage 2, *Nature Human Behaviour*).

The Stage 1 protocol, accepted in principle, is registered at
<https://doi.org/10.6084/m9.figshare.29504894>.

## Studies

**Study 1** — single-session experiment. Lonely adults recruited from CloudResearch Connect
were randomly assigned to one of six conditions: text–text, AI text–human voice,
AI voice–human text, voice–voice, a participant-choice condition allowing switching between
voice–voice and text–text mid-session, or a control condition spending the same time doing
nothing. All interactions lasted 15 minutes. Loneliness was measured immediately before and
after the session with an adapted state version of the UCLA loneliness scale. Final sample
after exclusions: N = 1,260.

**Study 2** — 10-day longitudinal experiment. Lonely adults were randomly assigned to
interact daily via text (text–text), via voice (voice–voice), or to a measurement-only
control. Participants completed one 15-minute session per day for ten consecutive days;
missing a day disqualified a participant from continuing. Loneliness and preference for
interacting with AI versus humans (0 = strongly prefer humans, 100 = strongly prefer AI)
were measured before and after each session. Final sample: 461 participants completed all
ten days (text = 153, voice = 137, control = 171).

## Repository contents

```
e1/
  e1.R                                     Study 1 analysis script
  data.csv                                 Study 1 survey data (4,743 session records; 1,260 after exclusions)
  sampled_conversations_for_manual_coding.csv
                                           200 sampled sessions with two independent human
                                           codes for emotional content (classification_k,
                                           classification_z) and the LLM code (label).
                                           Semicolon-delimited. Contains no conversation text.
e2/
  e2.R                                     Study 2 analysis script
  d1.csv … d10.csv                         Daily survey data, one file per study day
lit/
  Supplemental Materials Table S1 …xlsx    Lasgaard et al. normative table, used for the
                                           percentile comparison reported in the SI
LAB_LOG.md                                 Data-collection dates for both studies
```

[LAB_LOG.md](LAB_LOG.md) records the collection window for every study day, all of which fall
after in-principle acceptance and after the Stage 1 registration.

Participants are identified only by a pseudonymous `worker_id` (e.g. `P03283`). The
crosswalk to the original platform identifiers is not part of this repository.

Note that `e1/data.csv` and `e2/d*.csv` are Qualtrics exports: the two rows directly beneath
the header hold the full question wording and the Qualtrics import metadata rather than
participant data, and the scripts drop them during cleaning. In both studies
`lonely_before_1..3` are the three loneliness items before a session. The post-session items
are `lonely_after_1..3` in Study 2; in Study 1 they are `lonely_after_interac_1..3` for the
interaction conditions and `lonely_after_nothing_1..3` for the do-nothing control.
`preference_before_1` / `preference_after_1` (Study 2) are the AI-versus-human preference
rating, `social_anxiety_1..17` the Study 1 social anxiety scale, and `condition` the assigned
condition.

## Reproducing the analyses

Requires R (tested on 4.4.1). Both scripts install their own dependencies through `pacman`
on first run, so no manual package installation is needed. Run each script from inside its
own directory, because all paths are relative:

```sh
cd e1 && Rscript e1.R
cd e2 && Rscript e2.R
```

Each script prints its full results to standard output and writes figures to a `plots/`
subdirectory that it creates. Redirect the output to a file to browse it, for example
`Rscript e1.R > e1_output.txt`. Both scripts run start to finish from the data in this
repository and reproduce the sample sizes and test statistics reported in the manuscript,
except the conversation-based analyses listed under "Data not included" below.

Two switches control the Study 2 robustness analyses. `GET_ONLY_FINISHED` (near the top of
`e2.R`) defaults to `TRUE`, restricting the sample to the 461 participants who completed all
ten days, which is the basis of every number in the main text; setting it to `FALSE` and
re-running reproduces the "Full-sample replication" section of the SI. Re-running with the
flag flipped overwrites the figures in `e2/plots`, so restore them afterwards if you need the
manuscript versions. The matched-sample replication (H1–H4) and the propensity-score matching
appendix both run automatically on a default run; no flag is required.

## Data not included

Individual-level conversation data — the voice transcripts and text messages exchanged
between participants and the AI companion, and the underlying chatbot message logs and voice
recordings — are **not** included. Our IRB-approved protocol does not permit sharing this
conversational content, given its sensitive and potentially re-identifying nature.

The analysis scripts are written so that everything else still runs without those files. The
sections that depend on them detect their absence, print a `[reproducibility]` notice, and
skip, leaving the rest of the pipeline unaffected. Concretely, the following reported
analyses **cannot** be reproduced from this repository:

- Latent Dirichlet Allocation topic modelling of participants' messages (Study 1).
- Engagement metrics derived from the message logs — number of messages sent and average
  words per message — and the regression predicting loneliness reduction from them (Study 1).
- The participant word-count analyses added at a reviewer's request (not pre-registered):
  messages, words per message and total words sent by condition, and the models adjusting
  the modality comparisons for words sent (Study 1, SI Table A1; Study 2, SI Table B10).
- The modality-switching analyses for the participant-choice condition (Study 1).

Everything else reproduces, including all hypothesis tests in both studies, the mixed-effects
and latent growth models, the Bayes factors, the propensity-score matching appendix, and
Cohen's kappa for the two human coders of emotional content.
