# Build notes for the released PDF

The released `paper/main.pdf` is the byte output of one build of `paper/main.tex` with the two DOI
macros filled in:

```sh
cd paper
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

Tectonic `0.17.0` (Windows x86_64 release, `xtask -X compile`), TeX Live bundle fetched by Tectonic.
The fixed epoch makes the build reproducible: two independent builds in separate directories produced
identical bytes, SHA-256 `35b60cbbe38801a5a2e38f4dee8cf9b1ec025a7854bd11e1df770cf107613b57`,
259940 bytes, MD5 `82891f76893c631b358062c7c4da4977`, 24 pages.

Log statistics of the released build (`paper/main.log`, 16621 bytes):

- exit code `0`, `Output written on main.xdv (24 pages, 555060 bytes)`
- `0` overfull boxes
- `3` underfull boxes: badness 1946 in paragraph at line 233, badness 2126 in paragraph at line 278,
  badness 10000 in a table row paragraph at lines 54 to 57. They are benign short lines in the two
  floor tables and in one equation display
- `0` undefined references or citations, `0` LaTeX or package warnings, `0` missing glyphs
- empty bibliography log apart from the plain BibTeX banner

The released PDF carries exactly three external links, extracted from its own link annotations:
`https://doi.org/10.5281/zenodo.24428099`, `https://doi.org/10.5281/zenodo.24428101` and
`https://orcid.org/0000-0001-5839-5197`.

The source of the manuscript is the reviewed text admitted with SHA-256
`2ff2de5a756d618cf3fa064ad55f96dac39defb5d9e5739460fb773f085a226b`. The released `paper/main.tex`
differs from it on exactly two lines, the definitions of `\paperdoi` and `\suppdoi`. Its SHA-256 is
`aaf3cc86581b91c2294345d57700576e40017257ae1d2b976b95ff894a794dd4`. `paper/refs.bib` is unchanged,
SHA-256 `0f5a4df474e1e0189821028a452446b83973aab9f4f1cf15fef938e82e179c30`.
