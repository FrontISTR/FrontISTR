# Review Criteria for FrontISTR

This document sets out what to check when reviewing a Merge Request (MR), and how to proceed with the review.

## Basic stance

- Review is done to work out, together with the submitter, how a valuable change can be brought into FrontISTR.
- When requesting a change, give the reason unless it is obvious.
- Distinguish what must be addressed before merging from what can be improved later. Ensure the necessary quality, but do not ask for work that the MR does not need.
- Do not hold an individual reviewer responsible for oversights or remaining issues. Problems found after merging, and issues left open by an MR, are handled as new Issues.
- Whether to accept a feature that reviewers cannot decide on alone, or to remove a feature after it has been merged, is finally decided by the Software WG.

## What to check in review

- Check the MR from the aspects below, and request changes if there are problems that need to be addressed.
- Make clear what you have checked and what you could not judge, and share the checking with other reviewers as needed.
- Keep the scope of requested changes proportionate to the impact of the change. For the level of change requests, refer to the actual review examples.
- In a change request, show the specific problem and its impact, and distinguish problems introduced by this change from problems that already existed.

| Aspect | What to check | Examples from actual MRs and reviews |
|---|---|---|
| Correctness | Is the change consistent with the intended behavior and the mathematics? Does it also work correctly under related analysis conditions, in parallel execution, and in error handling? | [Pointed out that results change for the same physical model](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/700#note_4080960348) / [Pointed out that the coordinate system was not passed to the internal force and stress](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/699#note_4068021673) / [Pointed out that some ranks skip communication](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/576#note_3824847816) |
| Verification | Is the verification proportionate to the change? Do the added tests exercise the target feature, and can they detect bugs and unintended changes? | [Pointed out a test that passes even with the fix reverted](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/706#note_4080839899) / [Pointed out that the added tests were not run in CI](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/689#note_4081518373) / [Pointed out a compression test that does not go through the compression feature](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) |
| Design and division of responsibilities | Are the responsibilities and layers of each part kept, without details that belong to another layer creeping in? Is the change consistent with the existing patterns used in similar processing? | [Separated contact-specific processing from time increment control](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/691#note_4041722311) / [Separated assembly from element computation](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/593#note_3828414349) / [Pointed out analysis conditions brought into the middleware](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/642#note_3885639989) |
| Input/output specification | Does the change follow the existing item names and formats, and the existing distinctions in what each kind of item is used for? Are there problems with compatibility or with the impact on users? | [Pointed out specifying the hardware type in the analysis input](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/597#note_3751145248) / [Moved the damping coefficient from an option to a data line](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/577) |
| Complexity and extensibility | Does the change make the code harder to understand or to change later, through excessive branches, duplication, or dependencies, or through unclear grouping of processing? Does it add unnecessary burden to other features? | [Pointed out double management of stiffness matrix construction](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/504#note_2845420845) / [Asked to gather conditional branches in one place instead of scattering them per element](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/643#note_3827680479) / [Consolidated material processing and clarified responsibilities (a good example)](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/706#note_4080839899) |
| Consistency between description and implementation | Do comments, logs, the MR description, and the manual match the implementation? Are the scope of support and the impact on existing input explained? | [Pointed out descriptions and computed results left outdated after a fix](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/700#note_4094847272) / [Pointed out a mismatch between the log and the manual](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) |
| Intellectual property | For code, comments, test data, and so on, are there concerns involving third-party copyrights or patents, such as improper use of material from other software? This is checked when a reviewer notices a suspicion; an exhaustive investigation is not required. | [Pointed out a verbatim translation from CSparse and a missing attribution for LAPACK](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) / [Rechecked missing license notices](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3946065466) |

### Notes on design and division of responsibilities

- The basis for design is the structure and intent that can be read from the existing code and comments, and the patterns used in similar processing. It is not assumed that the design is written down in a design document.
- Details without a clear existing pattern or description are left to the judgment of the implementer and the reviewers.
- When a new design policy is decided, record it in the structure of the code or in comments, so that it can be read in later implementation and review.

## What is not a main target of review

- **Minor coding conventions, indentation, and the like**  
  This does not prevent pointing them out or requesting changes, but the aim is to unify them through automated checks and similar means.
- **Usefulness or priority of the feature itself**  
  For a feature that meets the criteria above, an individual reviewer's agreement is not a condition for merging. If someone doubts the need for a feature, the principle is that they open a separate Issue to discuss it. However, a reviewer may decline a change that is clearly outside the scope of FrontISTR.

## Changes left to the Software WG

When reviewers alone cannot decide whether to merge a change, or when a change needs a new policy decision on shared design, specifications, or compatibility, as in the cases below, organize the points at issue and leave the decision to the Software WG. Fixes and extensions that follow existing policies and patterns proceed through normal review.

- Changes that alter the existing policy on the division of responsibilities or the dependencies of common parts
- Changes that place significant constraints on adding or modifying other features, where it must be decided whether those constraints are acceptable
- Changes that alter the meaning or compatibility of existing input or output, and the introduction of new I/O formats that do not follow existing patterns

Even when the decision is left to the WG, continuing the checks and change requests along the review aspects above, and making the points at issue and the options concrete, is welcome.
