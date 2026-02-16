# Figure Generation Prompt Template (PreFerS v2)

Use this template to generate future figures with consistent PreFerS visual language.

## Copy-ready template

```markdown
Create one [FIGURE_TYPE] for [AUDIENCE] with **content-only** output unless otherwise stated.

### Communication goal
- Main message: [MAIN_MESSAGE]
- Decision/use context: [MEETING_OR_REPORT_CONTEXT]
- What viewers should understand in 5 seconds: [FIVE_SECOND_TAKEAWAY]

### Content blocks (exact text or near-exact)
- Block A: [TEXT]
- Block B: [TEXT]
- Block C: [TEXT]
- Block D: [TEXT]
- Block E: [TEXT]

### Structure and layout
- Layout pattern: [FLOW | ARGUMENT MAP | COMPARISON GRID | 2x2 TRADEOFF | TIMELINE]
- Reading order: [LEFT_TO_RIGHT | TOP_TO_BOTTOM | RADIAL]
- Group bands/sections: [SECTION_LABELS]
- Connector logic: [ARROW_RULES]
- Text density: [MAX_LINES_PER_BLOCK]

### Visual system (mandatory)
- Canonical source palette (original PreFerS):
  `#191538 #3C4C98 #2C80C4 #234966 #1B8A4D #8BC53F #DDE653 #F5CA0C #EDA211`
- Render variant for fills/cards: `PreFerS_softlight`
- Overall style: soft-light, antiqued, simplified Material-inspired
- Background: warm off-white with subtle paper texture
- Borders: thin and low-contrast
- Shadows: shallow and soft
- Connectors: dark desaturated blue-gray

### Legibility constraints
- High contrast for text at presentation scale
- Keep each text block concise ([2-4] lines recommended)
- Avoid clutter; prioritize hierarchy and spacing
- Icons optional; if used, minimal and non-dominant

### Output constraints
- Do [NOT_INCLUDE_ELEMENTS] (e.g., no title, no takeaway quote, no footnote)
- Target placement: [SLIDE_AREA] in 16:9 deck
- Visual intent: clear to general readers, still credible to technical audience
```

## Quick guide: choosing the best figure type
- **Argument Map**: best for “why we chose X” persuasion with multiple reasons.
- **Flowchart**: best for process sequence and causality.
- **Comparison Grid**: best for alternatives/options.
- **2x2 Tradeoff**: best for prioritization decisions.
- **Timeline**: best for phase plans and milestones.

## Recommended default
If unsure, use **Argument Map** with 1 center claim + 4–6 evidence cards + 1 synthesis strip.
