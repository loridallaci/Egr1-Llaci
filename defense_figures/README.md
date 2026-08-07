# defense_figures/

Figures used in the thesis defense deck **`LoridaLlaci_ThesisDefense_082526.pptx`**
(OneDrive copy, modified 2026-08-04 14:27 — this is the newer of the two copies;
a stale copy also exists at `Documents\Egr1\`), consolidated here on 2026-08-06 so
the analyses behind the talk live with the code that made them.

## How this was built

Every image embedded in the deck was MD5-hashed and matched against ~5,000
image/PDF files under the repo, `Final Submission/`, and `Egr1CC_vs_Egr1KDBulkRNA_FINAL_July2025/`.

- deck slides: **92**
- distinct embedded images: **132**
- exactly matched to a file on disk: **30**
- of those, already tracked in this repo: **7**
- copied into this folder: **23** (prefixed `s<NN>_` with the slide number)
- **untraced: 102**

## Files

| File | Contents |
|---|---|
| `manifest_matched.csv` | slide number, slide title, filename, whether it was already in the repo, its repo path, and the original location it was copied from |
| `manifest_untraced.csv` | the 102 deck images that could **not** be matched to a source file, with slide number and title |

## Why 102 could not be traced

PowerPoint re-encodes most pasted images (and rasterises vector paste-ins), so the
bytes in the deck no longer match the bytes of the file that produced them. An
unmatched entry does **not** mean the source file is missing — only that it could not
be identified automatically. The extracted images are available for manual matching at:

```
<scratchpad>\deck_unmatched_images\
```

named `s<NN>_imageNN.<ext>` so each can be compared against the slide it appears on.

## Related canonical outputs also brought into the repo

- `09_drug_screen/Fig5_revision/` — the 2026-07-28 Figure 5 panel set (B–E plus the
  Egr1-dependence pie). This is the **raw-p** version: panel B's hit counts are
  WT 61 male-biased + 2 female-biased = 63, KD 51 + 2 = 53.
  The older `09_drug_screen/Figure5_panelB_sexdiff_scatter.*` (2026-07-20) uses **BH q**
  and is superseded.
- `tables/SupplementaryTables_Final.xlsx` — added earlier on branch `add-supplementary-tables`.

## Known inconsistency to resolve

Two versions of the Egr1-dependence pie disagree on the denominator:

- `09_drug_screen/Egr1dependence_pie_n64.pdf` (08-04): 37 dependent + 27 independent = **n=64**
- `09_drug_screen/Fig5_revision/panelC_pie_Egr1dependence.pdf` (07-28): 37 + 26 = **n=63**

The difference is **THIOGUANINE**, which appears twice in the 85-compound plate
(85 rows, 84 unique names) and is plotted as "Thioguanine (run 1)" and
"Thioguanine (run 2)" in the dumbbell. n=64 counts both runs; n=63 collapses them.
Only one convention should reach the thesis or the paper.
