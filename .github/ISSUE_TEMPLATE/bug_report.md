---
name: Bug report
about: A wrong value, a missing or unexpected result, an error or a crash
labels: bug
---

**What happened, and what you expected**
A sentence or two. For a wrong value: which variable, where (time, position,
height) and why you think it is wrong (e.g. another instrument, a model,
physical limits).

**The call**
```python
# the loader call, exactly as run
```

**The output's attributes**
Paste the result of `print(ds.attrs)` (it records the call's options, the
package version, the files read, their data versions, files not found and
the data states). If the call raised an error, paste the whole traceback
instead.

```text

```

**Where it ran**
Operating system, Python version, and how the package was installed (uv,
pip, which version or commit).

**Anything else**
Whether it happens for other dates or instruments, or only for one file.
