# Toy Workflow Example

This directory provides a minimal structural example for the active UME
workflow. The files are intentionally tiny and are primarily meant to document
the expected file layout and naming conventions.

Because the full UME workflow depends on Perl and domain-specific input
formats, this toy example is best treated as a smoke-test scaffold rather than
a biologically meaningful run.

## Files

- `files.mergein.lst`
  Example source list used by `clist` and `discovery`
- `in_prod_mpileup_files.lst`
  Example production key file
- `toy-discovery.ll`
  Example discovery site-definition file

## Example commands

```bash
bin/ume clist \
  -f /path/to/All.tlone.Merge.db.gz \
  -k examples/toy-workflow/files.mergein.lst \
  -r /path/to/reference.fa.OL
```

```bash
bin/ume discovery \
  -f /path/to/All.tlone.Merge.db.gz \
  -k examples/toy-workflow/files.mergein.lst \
  -c 1 \
  -l 5 \
  -r /path/to/reference.fa.OL
```

```bash
bin/ume production \
  -f examples/toy-workflow/in_prod_mpileup_files.lst \
  -i /path/to/mpileup_dir \
  -k examples/toy-workflow/toy-discovery.ll \
  -r /path/to/reference.fa.OL \
  -d /tmp/ume-toy-output
```
