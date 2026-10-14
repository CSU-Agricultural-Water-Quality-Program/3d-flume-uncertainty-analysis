# Publishing a citable release

The repository is prepared for a first release, **v1.0.0**. A Zenodo DOI has
not been assigned. Do not insert an invented or example DOI into the manuscript.

Citation authors follow the manuscript's `ASABE_AuthorInformation_MD.docx`
in `1_First Submission/`: Deleon, Cabot, Brown, Ladow, White, Schlatter, and
Wardle, in that order. Names and affiliations were read from that form;
the manuscript and author-information document were not edited. No ORCID
identifiers are included because none were supplied in that form.

## Publish through the GitHub–Zenodo integration

1. Review `CITATION.cff` for author names/order, affiliation, version, and the
   existing GPL version 2 license. Add author-supplied ORCID URLs if desired.
   Set `date-released` to the actual publication date when publishing.
2. Commit and push the release preparation, including removal of tracked
   scratch files. `.gitignore` prevents future scratch files from being added;
   `.gitattributes` excludes them from source archives. The local files remain.
3. Sign in to [Zenodo](https://zenodo.org/), connect your GitHub account, and
   [enable this repository](https://help.zenodo.org/docs/github/enable-repository/).
   For the organization repository, the integration may need organization access.
4. Create a GitHub release at the reviewed commit, tagged **v1.0.0**, with the
   title **v1.0.0 — Metric flume calibration and uncertainty analysis**. Use
   `CHANGELOG.md` as the release description, including the diagnostic status.
   Publish a normal release rather than leaving it as a draft.
5. Wait for Zenodo to process the release, then check the archived files,
   authors, title, version, license, and DOI on the record. This follows
   [Zenodo's release workflow](https://help.zenodo.org/docs/github/archive-software/github-upload/).
6. Cite the **version DOI** in the paper so readers retrieve the exact code used.
   Add the assigned DOI to `CITATION.cff` as `doi: "10.5281/zenodo.…"` and link
   the record from the README. The ellipsis here is explanatory, not a DOI to copy.
   That metadata follow-up does not alter the already archived release files.

`CITATION.cff` is the single citation metadata file. GitHub exposes a citation
button and Zenodo can import these metadata. No `.zenodo.json` is needed for
this release; if one is added later, it overrides CFF metadata in Zenodo.
See [Zenodo's CFF guidance](https://help.zenodo.org/docs/github/describe-software/citation-file/).

## Files to archive

Include the main R script and every supporting module, optional tools, input
CSVs and provenance, reference documents, saved fit/diagnostics, historical
scripts/outputs, README files, citation metadata, changelog, and license.
Exclude `tmp/`, `release/`, `.Rhistory`, `.RData`, `.Ruserdata`, and `.Rproj.user/`.
The manuscript snapshots, extracted manual text, and compiled executables in
`tmp/` are local preparation artifacts, not part of the citable release.

After committing, an optional manual-upload archive can be created without
running R:

```powershell
New-Item -ItemType Directory -Force release | Out-Null
git archive --format=zip --prefix=3d-flume-uncertainty-analysis-1.0.0/ --output=release/3d-flume-uncertainty-analysis-1.0.0.zip HEAD
```

If uploading manually, select **Software**, use the CFF title/authors/version
and existing license, and upload this ZIP. Zenodo also lets you
[reserve a DOI before publication](https://help.zenodo.org/docs/deposit/describe-records/reserve-doi/).
Use either the GitHub integration or a manual deposit for a given release to
keep its citation unambiguous.

## Manuscript citation

Once the record exists, export the citation from Zenodo for your reference
manager. Its form will be:

> Deleon, E., Cabot, P., Brown, A. J., Ladow, J. F., White, D. C., Schlatter, K.,
> & Wardle, E. (2026). *3D-Printed Trapezoidal Flume Calibration and Bayesian
> Uncertainty Analysis* (Version 1.0.0) [Computer software]. Zenodo.
> [The version-specific DOI assigned to this release]

Suggested availability sentence, after publication:

> The calibration scripts and supporting analysis files are archived in Zenodo
> (version 1.0.0; insert the assigned version DOI).

Future code/data changes receive a new version and release. Zenodo links the
versions while providing a distinct identifier for each; see
[version management](https://help.zenodo.org/docs/deposit/manage-versions/).
