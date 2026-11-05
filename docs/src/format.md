
**`docs/src/format.md`**
```markdown
# File Format

`DiFloT.jl` uses the standard DIMACS `.min` format. However, instead of monetary costs, the last field on arc lines represents the **transit time** `τ_a`.

```text
c This is a comment
p min 3 3
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3