#!/usr/bin/env python3
"""
Load and parse Wind bow shock crossings from Excel database.
Provides easy access to crossing events, passes, and time windows.
"""

import openpyxl
import datetime as dt
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def _parse_datetime(value) -> Optional[dt.datetime]:
    """
    Parse datetime from various formats.
    Handles: datetime objects, ISO strings, other string formats.
    Returns None if value is None or empty.
    """
    if value is None:
        return None

    # Already a datetime object
    if isinstance(value, dt.datetime):
        if value.tzinfo is None:
            return value.replace(tzinfo=dt.timezone.utc)
        return value

    # String format
    if isinstance(value, str):
        s = str(value).strip()
        if not s:
            return None

        # Try ISO format first
        try:
            parsed = dt.datetime.fromisoformat(s)
            if parsed.tzinfo is None:
                parsed = parsed.replace(tzinfo=dt.timezone.utc)
            return parsed
        except ValueError:
            pass

        # Try common datetime formats
        for fmt in [
            "%Y-%m-%d %H:%M:%S.%f",
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%dT%H:%M:%S.%f",
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%d %H:%M",
            "%Y-%m-%dT%H:%M",
        ]:
            try:
                parsed = dt.datetime.strptime(s, fmt)
                if parsed.tzinfo is None:
                    parsed = parsed.replace(tzinfo=dt.timezone.utc)
                return parsed
            except ValueError:
                continue

        # If all else fails, raise error
        raise ValueError(f"Could not parse datetime: {value}")

    return None


class BowShockCrossing:
    """Represents a single bow shock crossing event."""

    def __init__(self, pass_id: int, sc_direction: str, time: dt.datetime,
                 crossing_direction: str, broad_start: Optional[dt.datetime],
                 broad_end: Optional[dt.datetime], zoom_start: Optional[dt.datetime],
                 zoom_end: Optional[dt.datetime]):
        """
        Initialize a crossing event.

        Args:
            pass_id: Pass/event ID
            sc_direction: Spacecraft location ('in' or 'out')
            time: Crossing timestamp
            crossing_direction: Crossing direction ('in' or 'out')
            broad_start: Broad window start time
            broad_end: Broad window end time
            zoom_start: Zoomed window start time
            zoom_end: Zoomed window end time
        """
        self.pass_id = pass_id
        self.sc_direction = sc_direction
        self.time = time
        self.crossing_direction = crossing_direction
        self.broad_start = broad_start
        self.broad_end = broad_end
        self.zoom_start = zoom_start
        self.zoom_end = zoom_end

    def has_broad_window(self) -> bool:
        """Check if this crossing has a broad time window."""
        return self.broad_start is not None and self.broad_end is not None

    def has_zoom_window(self) -> bool:
        """Check if this crossing has a zoomed time window."""
        return self.zoom_start is not None and self.zoom_end is not None

    def duration_broad(self) -> Optional[float]:
        """Get broad window duration in hours."""
        if self.has_broad_window():
            delta = self.broad_end - self.broad_start
            return delta.total_seconds() / 3600
        return None

    def duration_zoom(self) -> Optional[float]:
        """Get zoom window duration in minutes."""
        if self.has_zoom_window():
            delta = self.zoom_end - self.zoom_start
            return delta.total_seconds() / 60
        return None

    def __repr__(self):
        return (f"BowShockCrossing(pass={self.pass_id}, time={self.time}, "
                f"direction={self.crossing_direction}, sc_loc={self.sc_direction})")

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            'pass_id': self.pass_id,
            'sc_direction': self.sc_direction,
            'time': self.time,
            'crossing_direction': self.crossing_direction,
            'broad_start': self.broad_start,
            'broad_end': self.broad_end,
            'zoom_start': self.zoom_start,
            'zoom_end': self.zoom_end
        }


class BowShockDatabase:
    """Load and manage bow shock crossings database."""

    def __init__(self, xlsx_path: Optional[Path] = None):
        """
        Initialize database.

        Args:
            xlsx_path: Path to bs_crossings_V02.xlsx. If None, searches relative to this file.
        """
        if xlsx_path is None:
            # Try to find file relative to this module
            script_dir = Path(__file__).parent
            xlsx_path = script_dir / "bs_crossings_V02.xlsx"

            # Fallback to parent directory
            if not xlsx_path.exists():
                xlsx_path = script_dir.parent / "bs_crossings_V02.xlsx"

            # Repository layout: database_files/ beside this module
            if not xlsx_path.exists():
                xlsx_path = script_dir / "database_files" / "bs_crossings_V02.xlsx"

        self.xlsx_path = Path(xlsx_path)
        self.crossings: List[BowShockCrossing] = []
        self.passes: Dict[int, List[BowShockCrossing]] = {}

        self.load()

    def load(self) -> None:
        """Load crossing data from Excel file."""
        if not self.xlsx_path.exists():
            raise FileNotFoundError(f"Database not found: {self.xlsx_path}")

        wb = openpyxl.load_workbook(self.xlsx_path, data_only=True)
        ws = wb['bs_crossings_V02']

        self.crossings = []
        self.passes = {}

        for row in ws.iter_rows(min_row=2, values_only=True):
            # Stop at first empty row
            if row[0] is None:
                break

            # Parse row
            pass_id = int(row[0])
            crossing_direction = row[1]  # Column B: crossing direction (in/out)
            crossing_time = _parse_datetime(row[2])  # Column C: crossing timestamp
            sc_direction = row[3]  # Column D: spacecraft location (in/out)
            broad_start = _parse_datetime(row[4])  # Column E: broad window start
            broad_end = _parse_datetime(row[5])  # Column F: broad window end
            zoom_start = _parse_datetime(row[6])  # Column G: zoom window start
            zoom_end = _parse_datetime(row[7])  # Column H: zoom window end

            # Create crossing object
            crossing = BowShockCrossing(
                pass_id=pass_id,
                sc_direction=sc_direction,
                time=crossing_time,
                crossing_direction=crossing_direction,
                broad_start=broad_start,
                broad_end=broad_end,
                zoom_start=zoom_start,
                zoom_end=zoom_end
            )

            self.crossings.append(crossing)

            # Group by pass
            if pass_id not in self.passes:
                self.passes[pass_id] = []
            self.passes[pass_id].append(crossing)

        wb.close()

    def get_passes(self) -> List[int]:
        """Get sorted list of all pass IDs."""
        return sorted(self.passes.keys())

    def get_crossings_for_pass(self, pass_id: int) -> List[BowShockCrossing]:
        """Get all crossings for a specific pass."""
        return self.passes.get(pass_id, [])

    def get_crossings_with_windows(self) -> List[BowShockCrossing]:
        """Get only crossings that have at least one time window defined."""
        return [c for c in self.crossings if c.has_broad_window() or c.has_zoom_window()]

    def get_crossings_by_direction(self, direction: str) -> List[BowShockCrossing]:
        """Get crossings in specific direction ('in' or 'out')."""
        return [c for c in self.crossings if c.crossing_direction == direction]

    def get_closest_crossing(self, target_time: dt.datetime) -> BowShockCrossing:
        """Find crossing closest to target time."""
        return min(self.crossings,
                  key=lambda c: abs((c.time - target_time).total_seconds()))

    def summary(self) -> str:
        """Generate summary statistics."""
        total_crossings = len(self.crossings)
        total_passes = len(self.passes)
        with_windows = len(self.get_crossings_with_windows())
        inbound = len(self.get_crossings_by_direction('in'))
        outbound = len(self.get_crossings_by_direction('out'))

        min_time = min(c.time for c in self.crossings)
        max_time = max(c.time for c in self.crossings)

        return f"""
╔════════════════════════════════════════════════╗
║  Bow Shock Crossings Database Summary          ║
╚════════════════════════════════════════════════╝
Total Crossings:        {total_crossings}
Total Passes:           {total_passes}
Crossings w/ Windows:   {with_windows}
Inbound (entering):     {inbound}
Outbound (leaving):     {outbound}
Time Range:             {min_time} to {max_time}
"""

    def __repr__(self):
        return f"BowShockDatabase({len(self.crossings)} crossings, {len(self.passes)} passes)"

    def __len__(self):
        return len(self.crossings)


# Convenience function for notebook usage
def load_crossings(xlsx_path: Optional[Path] = None) -> BowShockDatabase:
    """
    Load bow shock crossings database (convenience wrapper).

    Usage in notebook:
        from bs_crossings_loader import load_crossings
        db = load_crossings()
        print(db.summary())
    """
    return BowShockDatabase(xlsx_path)


if __name__ == "__main__":
    # Test the loader
    print("Loading bow shock crossings database...")
    db = load_crossings()

    print(db.summary())
    print(f"\nDatabase: {db}")
    print(f"\nPasses: {db.get_passes()[:10]}... ({len(db.get_passes())} total)")

    # Show first pass crossings
    first_pass = db.get_passes()[0]
    crossings = db.get_crossings_for_pass(first_pass)
    print(f"\n{len(crossings)} crossings in Pass {first_pass}:")
    for i, c in enumerate(crossings[:5], 1):
        print(f"  {i}. {c}")
