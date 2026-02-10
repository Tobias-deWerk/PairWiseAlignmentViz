from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Union

AUTO_WIDTH_TOKEN = "auto"
DEFAULT_ANNOTATION_THICKNESS = 10.0
DEFAULT_ANNOTATION_ALPHA = 0.85
DEFAULT_REF_ANNOTATION_COLOR = "#984ea3"
DEFAULT_QUERY_ANNOTATION_COLOR = "#ff7f00"
DEFAULT_ANNOTATION_MAX_LAYERS = 3
DEFAULT_ANNOTATION_SPACING = 0.8
DEFAULT_ANNOTATION_LABEL_JITTER = 20.0


WidthArg = Union[float, str]


@dataclass(frozen=True)
class RenderParams:
    width: WidthArg = AUTO_WIDTH_TOKEN
    height: float = 6.0
    dpi: int = 150
    min_gap_size: int = 10
    block_size: int = 100
    min_sequence_identity: float = 0.7
    window_size: int = 20
    tick_interval: int = 10_000
    backbone_gap: float = 0.5
    backbone_thickness: float = 3.0
    bump_scale: float = 0.3
    mismatch_line_width: float = 0.4
    gap_max_height: float = 2.0
    gap_width: float = 50.0
    gap_height_scale: float = 0.001
    indel_height_scale: float = 0.01
    gap_label_size: Optional[float] = None
    annotation_label_size: Optional[float] = 10.0
    annotation_thickness: float = DEFAULT_ANNOTATION_THICKNESS
    annotation_alpha: float = DEFAULT_ANNOTATION_ALPHA
    reference_annotation_color: str = DEFAULT_REF_ANNOTATION_COLOR
    query_annotation_color: str = DEFAULT_QUERY_ANNOTATION_COLOR
    annotation_label_jitter: float = DEFAULT_ANNOTATION_LABEL_JITTER
    annotation_max_layers: int = DEFAULT_ANNOTATION_MAX_LAYERS
    annotation_spacing: float = DEFAULT_ANNOTATION_SPACING

    @classmethod
    def from_cli_args(cls, args: Any) -> "RenderParams":
        return cls(
            width=args.width,
            height=float(args.height),
            dpi=int(args.dpi),
            min_gap_size=int(args.min_gap_size),
            block_size=int(args.block_size),
            min_sequence_identity=float(args.min_sequence_identity),
            window_size=int(args.window_size),
            tick_interval=int(args.tick_interval),
            backbone_gap=float(args.backbone_gap),
            backbone_thickness=float(args.backbone_thickness),
            bump_scale=float(args.bump_scale),
            mismatch_line_width=float(args.mismatch_line_width),
            gap_max_height=float(args.gap_max_height),
            gap_width=float(args.gap_width),
            gap_height_scale=float(args.gap_height_scale),
            indel_height_scale=float(args.indel_height_scale),
            gap_label_size=args.gap_label_size,
            annotation_label_size=args.annotation_label_size,
            annotation_thickness=float(args.annotation_thickness),
            annotation_alpha=float(args.annotation_alpha),
            reference_annotation_color=str(args.reference_annotation_color),
            query_annotation_color=str(args.query_annotation_color),
            annotation_label_jitter=float(args.annotation_label_jitter),
            annotation_max_layers=int(args.annotation_max_layers),
            annotation_spacing=float(args.annotation_spacing),
        )

    @classmethod
    def from_payload(cls, payload: Dict[str, Any]) -> "RenderParams":
        def require(name: str, default: Any) -> Any:
            return payload.get(name, default)

        return cls(
            width=parse_width(require("width", AUTO_WIDTH_TOKEN)),
            height=to_float(require("height", 6.0), positive=True, name="height"),
            dpi=to_int(require("dpi", 150), positive=True, name="dpi"),
            min_gap_size=to_int(require("min_gap_size", 10), positive=True, name="min_gap_size"),
            block_size=to_int(require("block_size", 100), positive=True, name="block_size"),
            min_sequence_identity=to_float(require("min_sequence_identity", 0.7), min_value=0.0, max_value=1.0, name="min_sequence_identity"),
            window_size=to_int(require("window_size", 20), positive=True, name="window_size"),
            tick_interval=to_int(require("tick_interval", 10_000), min_value=0, name="tick_interval"),
            backbone_gap=to_float(require("backbone_gap", 0.5), min_value=0.0, name="backbone_gap"),
            backbone_thickness=to_float(require("backbone_thickness", 3.0), min_value=0.0, name="backbone_thickness"),
            bump_scale=to_float(require("bump_scale", 0.3), min_value=0.0, name="bump_scale"),
            mismatch_line_width=to_float(require("mismatch_line_width", 0.4), min_value=0.0, name="mismatch_line_width"),
            gap_max_height=to_float(require("gap_max_height", 2.0), min_value=0.0, name="gap_max_height"),
            gap_width=to_float(require("gap_width", 50.0), min_value=0.0, name="gap_width"),
            gap_height_scale=to_float(require("gap_height_scale", 0.001), min_value=0.0, name="gap_height_scale"),
            indel_height_scale=to_float(require("indel_height_scale", 0.01), min_value=0.0, name="indel_height_scale"),
            gap_label_size=parse_optional_size(require("gap_label_size", "NULL"), "gap_label_size"),
            annotation_label_size=parse_optional_size(require("annotation_label_size", 10.0), "annotation_label_size"),
            annotation_thickness=to_float(require("annotation_thickness", DEFAULT_ANNOTATION_THICKNESS), min_value=0.0, name="annotation_thickness"),
            annotation_alpha=to_float(require("annotation_alpha", DEFAULT_ANNOTATION_ALPHA), min_value=0.0, max_value=1.0, name="annotation_alpha"),
            reference_annotation_color=str(require("reference_annotation_color", DEFAULT_REF_ANNOTATION_COLOR)),
            query_annotation_color=str(require("query_annotation_color", DEFAULT_QUERY_ANNOTATION_COLOR)),
            annotation_label_jitter=to_float(require("annotation_label_jitter", DEFAULT_ANNOTATION_LABEL_JITTER), min_value=0.0, name="annotation_label_jitter"),
            annotation_max_layers=to_int(require("annotation_max_layers", DEFAULT_ANNOTATION_MAX_LAYERS), positive=True, name="annotation_max_layers"),
            annotation_spacing=to_float(require("annotation_spacing", DEFAULT_ANNOTATION_SPACING), min_value=0.0, name="annotation_spacing"),
        )


def parse_width(value: Any) -> WidthArg:
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == AUTO_WIDTH_TOKEN:
            return AUTO_WIDTH_TOKEN
    parsed = to_float(value, positive=True, name="width")
    return parsed


def parse_optional_size(value: Any, name: str) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, str) and value.strip().lower() in {"na", "null", ""}:
        return None
    return to_float(value, positive=True, name=name)


def to_float(
    value: Any,
    *,
    name: str,
    positive: bool = False,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a floating-point number") from exc
    if positive and parsed <= 0:
        raise ValueError(f"{name} must be positive")
    if min_value is not None and parsed < min_value:
        raise ValueError(f"{name} must be >= {min_value}")
    if max_value is not None and parsed > max_value:
        raise ValueError(f"{name} must be <= {max_value}")
    return parsed


def to_int(
    value: Any,
    *,
    name: str,
    positive: bool = False,
    min_value: Optional[int] = None,
) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if positive and parsed <= 0:
        raise ValueError(f"{name} must be positive")
    if min_value is not None and parsed < min_value:
        raise ValueError(f"{name} must be >= {min_value}")
    return parsed
