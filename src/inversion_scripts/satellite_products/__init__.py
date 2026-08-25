"""Satellite-product strategies used by IMI workflows."""

from .base import ObservationRequest, SatelliteProduct, SuperobservationResult
from .registry import get_satellite_product, supported_satellite_products

__all__ = [
    "ObservationRequest",
    "SatelliteProduct",
    "SuperobservationResult",
    "get_satellite_product",
    "supported_satellite_products",
]
