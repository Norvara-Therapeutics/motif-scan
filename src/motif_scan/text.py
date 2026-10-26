"""Number and sequence formatting shared by the terminal output and the report.

Only human-facing output is formatted here. The TSV and JSON outputs carry
unrounded values and full matched text, so a figure another program reads is
never the rounded one.
"""

from __future__ import annotations


def thousands(value: int) -> str:
    """Group an integer in thousands."""
    return f"{value:,}"


def plural(count: int, singular: str, suffix: str = "s") -> str:
    """Count something, agreeing the noun with the number."""
    return f"{thousands(count)} {singular}{'' if count == 1 else suffix}"


def number(value: float | None, places: int = 2, missing: str = "n/a") -> str:
    """Format an optional number to a fixed number of decimal places."""
    if value is None:
        return missing
    return f"{value:,.{places}f}"


def ratio(value: float | None, missing: str = "n/a") -> str:
    """Format an observed-to-expected ratio.

    Shown to two decimal places, and to one significant figure below 0.01,
    where two decimal places would print every rare motif as 0.00 and so lose
    the only part of the figure that carries information.
    """
    if value is None:
        return missing
    if 0 < value < 0.01:
        return f"{value:.0e}"
    return f"{value:,.2f}"


def residues(count: int, unit: str) -> str:
    """Count residues in the unit of their alphabet."""
    return f"{thousands(count)} {unit}"


def span(start: int, end: int) -> str:
    """Write a 1-based inclusive range the way positions are reported."""
    return f"{thousands(start)}..{thousands(end)}"


def truncate(text: str, width: int) -> str:
    """Shorten text to a width, marking that it was shortened.

    Used for a pattern or a description in a fixed-width table. The full text
    is always in the TSV and JSON outputs, so nothing is lost, only abbreviated.
    """
    if width <= 1 or len(text) <= width:
        return text
    return text[: width - 1] + "\u2026"


def wrap_sequence(sequence: str, width: int = 60) -> list[str]:
    """Break a sequence into fixed-width lines."""
    return [sequence[index: index + width] for index in range(0, len(sequence), width)]
