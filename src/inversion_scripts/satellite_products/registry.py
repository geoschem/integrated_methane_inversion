"""Registry of satellite products supported by the IMI."""

from __future__ import annotations

from src.inversion_scripts.satellite_products.base import SatelliteProduct
from src.inversion_scripts.satellite_products.products import (
    BlendedTropomiProduct,
    MethaneSatProduct,
    TropomiProduct,
)


_SATELLITE_PRODUCTS: dict[str, SatelliteProduct] = {
    product.name: product
    for product in (
        TropomiProduct(),
        BlendedTropomiProduct(),
        MethaneSatProduct(),
    )
}


def supported_satellite_products() -> tuple[str, ...]:
    """Return the canonical names of all registered products."""
    return tuple(sorted(_SATELLITE_PRODUCTS))


def get_satellite_product(name: str) -> SatelliteProduct:
    """Return the registered product named *name* or raise a useful error."""
    try:
        return _SATELLITE_PRODUCTS[name]
    except KeyError:
        supported = ", ".join(supported_satellite_products())
        raise ValueError(
            f"Unsupported satellite product {name!r}; supported products: {supported}"
        ) from None
