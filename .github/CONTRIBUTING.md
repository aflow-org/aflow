# Contributing Guidelines

### Contents

- [Asking Questions](#asking-questions)
- [Opening an Issue](#opening-an-issue)
- [Feature Requests](#feature-requests)
- [Submitting Pull Requests](#submitting-pull-requests)
- [Coding Style](#coding-style)

## Asking Questions

Questions should be asked in the [discussion](https://docs.github.com/en/discussions/quickstart) section of this repository. 
Please keep in mind that we don't provide formal support, and will not help with solving research questions. Questions 
should only focus on the usage of `aflow` itself.

## Opening an Issue
Before [creating an issue](https://help.github.com/en/github/managing-your-work-on-github/creating-an-issue), check if 
you are using the latest version of `aflow`. If you are not up-to-date, see if updating fixes your issue first.

### Reporting Security Issues
If you discover a security issue, please bring it to our attention right away! **DO NOT** file a public issue for 
security vulnerabilities. Instead, email us at **aflow_sec@materials.duke.edu**.

### Bug Reports and Other Issues
A great way to contribute to `aflow` is to send us a detailed issue when you encounter a problem. We always appreciate 
a well-written, thorough bug report. 

- **Review the [documentation]()** before opening a new issue.
- **Do not open a duplicate issue!** Search through existing issues to see if your issue has previously been reported. 
  If your issue exists, comment with any additional information you have. You may simply note "I have this problem too", 
  which helps prioritize the most common problems and requests.
- **Prefer using [reactions](https://github.blog/2016-03-10-add-reactions-to-pull-requests-issues-and-comments/)**, 
  not comments, if you simply want to "+1" an existing issue.
- **Be clear, concise, and descriptive.** Provide as much information as you can, including steps to reproduce, stack traces, 
  compiler errors, library versions, OS versions, and screenshots (if applicable).
- **Use [GitHub-flavored Markdown](https://help.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).** 
  Especially put code blocks and console outputs in backticks (```). This improves readability.

## Feature Requests
Feature requests are welcome! While we will consider all requests, we cannot guarantee your request will be accepted. 
We want to avoid [feature creep](https://en.wikipedia.org/wiki/Feature_creep). Your idea may be great, but also out-of-scope
for `aflow`. If accepted, we cannot make any commitments regarding the timeline for implementation and release. However, 
you are welcome to submit a pull request to help!

- **Do not open a duplicate feature request.** Search for existing feature requests first. If you find your feature 
  (or one very similar) previously requested, comment on that issue.
- Be precise about the proposed outcome of the feature and how it relates to existing features. Include implementation 
  details if possible.

## Submitting Pull Requests
While `aflow` is open-source it is not open-contribution. So all code changes will go through the maintainers hands before
being included in a new public release of `aflow`. Nevertheless, we **love** to receive pull requests! If we include
parts of a pull request we will add you as a contributor in the combined release commits, and your contribution will be
marked in at the relevant code sections.

Before [forking the repo](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) 
and [creating a pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) 
for non-trivial changes, it is usually best to first open an issue to discuss the changes, or discuss your intended 
approach for solving the problem in the comments for an existing issue.

*Note: All contributions will be licensed under [GPLv3](https://choosealicense.com/licenses/gpl-3.0/).*

- **Smaller is better.** Submit **one** pull request per bug fix or feature. A pull request should contain isolated 
  changes pertaining to a single bug fix or feature implementation. **Do not** refactor or reformat code that is 
  unrelated to your change. It is better to **submit many small pull requests** rather than a single large one. 
  Enormous pull requests will take enormous amounts of time to review, or may be rejected altogether.
- **Coordinate bigger changes.** For large and non-trivial changes, open an issue to discuss a strategy with the 
  maintainers. Otherwise, you risk doing a lot of work for nothing!
- **Prioritize understanding over cleverness.** Write code clearly and concisely. Remember that source code 
  usually gets written once and read often. Ensure the code is clear to the reader. The purpose and logic should be 
  obvious to a reasonably skilled developer, otherwise you should add a comment that explains it.
- **Follow existing coding style and conventions.** Keep your code consistent with the style, formatting, and 
  conventions in the rest of the code base. When possible, these will be enforced with a linter. Consistency makes 
  it easier to review and modify in the future.
- **Include test coverage.** Add unit tests when possible. Follow existing patterns for implementing tests.
- **Update the example project** if one exists to exercise any new functionality you have added.
- **Add documentation.** Document your changes with code doxy comments as described in [docs/README.md](/docs/README.md).
- **Use the repo's default branch.** Branch from and [submit your pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork) 
  to the repo's release branch.
- **[Resolve any merge conflicts](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-on-github)** that occur.
- When writing comments, use properly constructed sentences, including punctuation.
- Use spaces, not tabs.

## Coding Style

C++ code is formatted according to the `.clang-format` at the root of the repository. Your IDE will typically detect this file 
and apply formatting as you code, but you should also make sure you apply the formatting as you work. Your IDE will likely have
a keybind for formatting code, but you can also apply formatting with the `clang-format` program. If you have bash and awk 
available, you can use the `scripts/format-changed.sh` script to be notified of and apply formatting on changed lines since the 
base branch.

You should also use clang-tidy to check for warnings. You may tell your IDE to use the `.clang-tidy` or `.clang-tidy-lint` files 
at the root of the repository. The latter has less checks, but these are the ones we are stricter about. The checks in the `.clang-tidy-lint` 
should all pass before something can be merged. If you have bash and awk available, you can use the `scripts/tidy-changed.sh` script to 
be notified of and apply clang-tidy fixes on changed lines since the base branch.

Hint: Use `git diff --name-only` or `git diff --name-only --staged` to get a list of changed files to pass to run-clang-tidy and clang-format 
for whole-file changes, or use `scripts/format-changed.sh` and `scripts/tidy-changed.sh` to detect and apply fixes to lines changed since 
the base branch.

For reference, `scripts/format-changed.sh` and `scripts/tidy-changed.sh` will effectively be run before a merge can be completed. If your 
local linting messages/changes don't match those of the linter on the remote repository, check that you are using the same version of the 
llvm tools.

