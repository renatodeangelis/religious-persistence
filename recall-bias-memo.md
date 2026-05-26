---
title: "Add Health Recall Bias Validation — Tables"
date: "May 2026"
output: pdf_document
geometry: margin=1in
fontsize: 11pt
---

## Data and Measures

**Sample.** Wave 1 in-home interview ($N = 6{,}504$) linked to Wave 3 in-home
interview via respondent ID (AID). Inner join on AID yields $N = 4{,}882$
respondents present at both waves; 1,622 Wave 1 respondents are lost to
attrition by Wave 3. All Wave 3 respondents matched exactly to Wave 1.

**Variables.**

| Variable | Wave | Label | Role |
|---|---|---|---|
| PA22 | W1 | Parent self-reported religion | Ground truth for childhood religious environment |
| H1RE1 | W1 | Respondent current religion | Adolescent religious identity (~age 16) |
| H3RE1 | W3 | Respondent current religion | Adult religious identity (~age 22--26) |
| H3RE26 | W3 | Religion raised in | Retrospective recall of childhood religion |

**Category scheme.** A common 5-category scheme is constructed for comparability
across waves.

| Code | Category | Wave 1 mapping | Wave 3 mapping |
|---|---|---|---|
| 0 | None | Code 0 | Code 0 |
| 1 | Protestant | Codes 1--19, 27 (incl. JW, LDS, Christian Science, Unitarian) | Codes 1, 8 |
| 2 | Catholic | Code 22 | Code 2 |
| 3 | Jewish | Code 26 | Code 3 |
| 4 | Other | Codes 20, 21, 23--25, 28 | Codes 4--7 |

Refused, don't know, and not applicable codes are set to NA. Buddhist, Hindu,
and Muslim respondents are collapsed into Other due to thin cell counts.

---

**Table 1.** Child's own religion at Wave 1 (H1RE1) by PA22 missingness.

| Child religion (W1) | PA22 observed (*n*) | PA22 observed (%) | PA22 missing (*n*) | PA22 missing (%) |
|---|---|---|---|---|
| None | 456 | 10.7 | 84 | 13.7 |
| Protestant | 2,598 | 60.9 | 367 | 59.7 |
| Catholic | 939 | 22.0 | 109 | 17.7 |
| Jewish | 42 | 1.0 | 4 | 0.7 |
| Other | 153 | 3.6 | 41 | 6.7 |
| Missing | 79 | 1.9 | 10 | 1.6 |

---

**Table 2.** Child--parent religious discordance at Wave 1. Complete cases ($N = 4{,}188$).

| Parental religion (PA22) | *n* | Discordance rate |
|---|---|---|
| None | 265 | 0.434 |
| Protestant | 2,709 | 0.137 |
| Catholic | 1,042 | 0.211 |
| Jewish | 48 | 0.208 |
| Other | 124 | 0.573 |
| **Overall** | **4,188** | **0.188** |

---

**Table 3.** Child--parent religion cross-tabulation at Wave 1 (raw counts; rows = PA22, cols = H1RE1).

| | **None** | **Protestant** | **Catholic** | **Jewish** | **Other** |
|---|---|---|---|---|---|
| **None** | 150 | 84 | 19 | 1 | 11 |
| **Protestant** | 215 | 2,339 | 84 | 1 | 70 |
| **Catholic** | 70 | 131 | 822 | 2 | 17 |
| **Jewish** | 2 | 4 | 2 | 38 | 2 |
| **Other** | 19 | 40 | 12 | 0 | 53 |

---

**Table 4.** Retrospective recall (H3RE26) by actual parental religion (PA22). Row proportions; diagonal in **bold**. Overall mismatch rate = 0.253. Complete cases ($N = 4{,}209$).

| | **Recalled: None** | **Recalled: Protestant** | **Recalled: Catholic** | **Recalled: Jewish** | **Recalled: Other** | *Mismatch rate* |
|---|---|---|---|---|---|---|
| **PA22: None** | **0.508** | 0.265 | 0.125 | 0.004 | 0.098 | 0.492 |
| **PA22: Protestant** | 0.087 | **0.742** | 0.066 | 0.002 | 0.103 | 0.258 |
| **PA22: Catholic** | 0.062 | 0.053 | **0.877** | 0.000 | 0.009 | 0.123 |
| **PA22: Jewish** | 0.041 | 0.082 | 0.041 | **0.816** | 0.020 | 0.184 |
| **PA22: Other** | 0.137 | 0.476 | 0.137 | 0.008 | **0.242** | 0.758 |

---

**Table 5.** Mismatch type by adult religion at Wave 3 (H3RE1). Row percentages. Complete cases ($N = 3{,}812$).

| Adult religion (W3) | Concordant | More secular | More religious | Different religion |
|---|---|---|---|---|
| None | 52.3 | **31.5** | 4.3 | 11.8 |
| Protestant | 87.7 | 1.7 | 2.5 | 8.1 |
| Catholic | 83.6 | 1.2 | 2.2 | 12.9 |
| Jewish | 86.1 | 2.8 | 2.8 | 8.3 |
| Other | 31.4 | 4.1 | 4.6 | 60.0 |

---

**Table 6.** Logistic regression: P(recalled more secular) ~ adult religion + parental religion (+ adolescent religion). Outcome: recalled None when PA22 $\neq$ None. Protestant reference throughout. $N = 3{,}836$.

| | Model 1 | | Model 2 | |
|---|---|---|---|---|
| **Term** | **OR** | **95% CI** | **OR** | **95% CI** |
| *Adult religion (ref: Protestant)* | | | | |
| None | 38.85\*\*\* | [26.41, 57.15] | 31.28\*\*\* | [21.04, 46.50] |
| Catholic | 0.79 | [0.36, 1.73] | 1.17 | [0.53, 2.62] |
| Jewish | 5.39 | [0.49, 58.93] | 4.51 | [0.36, 56.01] |
| Other | 2.27\* | [1.18, 4.37] | 2.23\* | [1.16, 4.31] |
| *Parental religion (ref: Protestant)* | | | | |
| Catholic | 0.73 | [0.51, 1.04] | 1.15 | [0.73, 1.82] |
| Jewish | 0.24 | [0.04, 1.37] | 0.21 | [0.02, 2.27] |
| Other | 1.46 | [0.76, 2.81] | 1.21 | [0.58, 2.49] |
| *Adolescent religion at W1 (ref: Protestant)* | | | | |
| None | | | 3.25\*\*\* | [2.27, 4.65] |
| Catholic | | | 0.30\*\*\* | [0.17, 0.56] |
| Jewish | | | 1.60 | [0.14, 18.47] |
| Other | | | 1.21 | [0.61, 2.38] |

\* $p < .05$; \*\*\* $p < .001$
