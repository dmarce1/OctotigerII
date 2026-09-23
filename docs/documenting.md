# Documenting and formatting code {#documenting_code}

Use `///` for a short API description and `/** ... */` for equations, invariants,
parameters, or a longer explanation. Document what a caller must know: units,
array order, ownership, preconditions, failure behavior, and which synchronization
operation makes the data available. Describe a numerical method where it is
implemented, including deviations from its cited reference.

A linked author–year citation is ordinary Doxygen syntax:

```cpp
/** Two-wave HLL flux.
 * See @ref ref_harten1983 "Harten et al. (1983)".
 * @param normal Coordinate direction: 0=x, 1=y, 2=z.
 */
```

Add the complete record to `BIBLIOGRAPHY.md` using a heading anchor such as
`{#ref_harten1983}`. Follow the author–year and a/b/c rules at the top of that
file. `docs.sh` checks citation keys before running Doxygen. Doxygen warnings
fail the documentation build; inspect `warnings.log` in the chosen output folder.

Doxygen extracts the full API, including private implementation details, so a
reader can navigate all symbols. Detailed comments concentrate on contracts and
methods rather than restating each trivial operator. Extraction does not certify
that every symbol has prose documentation.

Keep one empty line between methods, functions, and class sections, and two
empty lines between larger declarations such as classes. `.clang-format`
preserves up to two empty lines but does not insert every required separator.
Use clang-format 21.1.8 with four-column tabs and K&R braces.

Use `class`, `typename` for template type parameters, and `using` aliases. Within
an access section, place constructors/destructors, member operators, member
functions, static functions, friend functions, and friend operators in that order
where declaration dependencies permit. Put `using std::sqrt;` and other needed
math declarations at the start of the function, then call the function unqualified
so argument-dependent lookup can find overloads.

Use `[[nodiscard]]` only when ignoring the result causes a specific correctness
hazard. An ignored calculation or getter that merely does nothing is not enough.
