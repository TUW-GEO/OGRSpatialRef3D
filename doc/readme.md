### Incompatibilities with PROJ.4 ###
This section describes the incompatibilities from existing coordinate transformation workflow (involving `OGRSpatialReference` and `OGRCreateCoordinateTransformation`).
The coordinate transformation in SpatialRef3D uses modified PROJ.4 `pj_transform` procedure. The algorithm which uses vertical grid shift from/to source/target SRS are replaced with height correction transform in `OGRSpatialReference3D`.

---

