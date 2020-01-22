# Examples

## As0

The `As0.output.json` file was created with the following command:

```bash
repp make seq -i As0.input.fa -o As0.output.json --addgene --settings twist.yaml -v
```

Which approximates to asking REPP to:

- make a plasmid sequence that's in `As0.input.fa`
- output build instructions in JSON to `As0.output.json`
- use source fragments from Addgene
- use synthesis settings from `twist.yaml`
- be verbose during the build
