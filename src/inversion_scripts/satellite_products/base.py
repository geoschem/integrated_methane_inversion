"""Common interfaces for satellite-product-specific behavior."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ObservationRequest:
    """Inputs that a satellite product may need to process one raw file."""

    filename: str | Path
    species: str
    start_date: Any
    end_date: Any
    xlim: tuple[float, float] | list[float]
    ylim: tuple[float, float] | list[float]
    use_water_observations: bool = False
    state_vector_path: str | Path | None = None
    gc_cache: str | Path | None = None
    time_threshold: str | None = None
    config: dict[str, Any] | None = None


@dataclass(frozen=True)
class SuperobservationResult:
    """Canonical superobservations and the shape of their model grid."""

    observations: np.ndarray
    grid_shape: tuple[int, ...]


class SatelliteProduct(ABC):
    """Strategy interface implemented by every supported retrieval product."""

    name: str
    goopy_raw_product_name: str
    visualization_source: str = "raw"

    @abstractmethod
    def observation_date(self, path: str | Path) -> datetime:
        """Return the acquisition date encoded in a raw or derived filename."""

    @abstractmethod
    def validate_directory(self, directory: Path) -> None:
        """Validate raw observation files in *directory*."""

    @abstractmethod
    def read_and_filter(
        self, request: ObservationRequest
    ) -> tuple[dict, np.ndarray] | None:
        """Read a raw file and return it with indices of valid observations."""

    @abstractmethod
    def create_superobservations(
        self, request: ObservationRequest
    ) -> SuperobservationResult | None:
        """Create canonical IMI superobservations for one raw file."""

    @abstractmethod
    def preview(self, request: ObservationRequest) -> dict | None:
        """Return product-neutral fields used by the IMI preview.

        The returned ``time`` field must contain datetime-like values rather
        than product-specific timestamp strings.
        """

    def download(self, start_date: str, end_date: str, destination: Path) -> None:
        """Download observations, or report that this product needs user data."""
        raise ValueError(
            f"{self.name} must be supplied with DataPathObs; "
            "automatic download is not supported."
        )
