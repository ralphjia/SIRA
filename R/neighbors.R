# =============================================================================
# R/neighbors.R
# Voxel geometry helpers and neighbor structure construction.
#
# Supports two data layouts:
#   Grid layout   : voxels on a d1 x d2 x d3 rectangular grid (env$coords NULL)
#   Coords layout : arbitrary spatial locations given as a V x 3 coordinate
#                   matrix (env$coords non-NULL, e.g. brain-masked MNI data)
#
# In both cases the same downstream structures are produced and stored in env:
#   env$neighbor_list    - list of length V; element v is integer vector of
#                          6-connected neighbor voxel indices for voxel v
#   env$edge_list        - 2-column integer matrix; undirected edges (i < j)
#   env$edge_structure   - list(edges, index) for fast SAD look-ups
#   env$smoothing_matrix - sparse V x V matrix used in initialization
# =============================================================================


# -----------------------------------------------------------------------------
# .sira_build_neighbors
# Dispatcher: calls the grid or coordinate variant based on env$coords.
# -----------------------------------------------------------------------------

#' @keywords internal
.sira_build_neighbors <- function(env) {
  if (!is.null(env$coords)) {
    .sira_build_neighbors_from_coords(env)
  } else {
    .sira_build_neighbors_from_grid(env)
  }
}


# -----------------------------------------------------------------------------
# .sira_build_neighbors_from_grid
# Original grid-based implementation (column-major voxel numbering).
# -----------------------------------------------------------------------------

#' @keywords internal
.sira_build_neighbors_from_grid <- function(env) {
  V  <- env$V
  d1 <- env$d1;  d2 <- env$d2;  d3 <- env$d3

  neighbor_list  <- vector("list", V)
  edge_list_rows <- vector("list", V)

  for (v in seq_len(V)) {
    nb <- .voxel_neighbors(v, d1, d2, d3)
    neighbor_list[[v]] <- nb
    new_edges <- nb[nb > v]
    if (length(new_edges) > 0L)
      edge_list_rows[[v]] <- cbind(v, new_edges)
  }

  .sira_finalise_neighbor_structures(neighbor_list, edge_list_rows, V, env)
}


# -----------------------------------------------------------------------------
# .sira_build_neighbors_from_coords
# Coordinate-based neighbor construction for irregular/masked data.
#
# Two voxels are 6-connected if they differ by exactly one voxel spacing in
# exactly one axis and are coincident in the other two.  The per-axis spacing
# is detected as the minimum positive gap between unique coordinate values.
# Coordinates are normalized to an integer grid before hashing, which handles
# floating-point representations of regular grids (e.g. MNI mm coordinates).
# -----------------------------------------------------------------------------

#' @keywords internal
.sira_build_neighbors_from_coords <- function(env) {
  V      <- env$V
  coords <- env$coords   # V x 3

  # --- Per-axis spacing detection ----------------------------------------
  spacing <- apply(coords, 2L, function(ax) {
    u <- sort(unique(ax))
    if (length(u) < 2L) return(NA_real_)
    g <- min(diff(u))
    if (g < 1e-6) NA_real_ else g
  })
  # Degenerate axes (all voxels share one coordinate) get a large sentinel
  # so offsets along that axis never match any real voxel.
  spacing[is.na(spacing)] <- .Machine$double.xmax

  # --- Normalize to integer grid -----------------------------------------
  grid <- round(sweep(coords, 2L, spacing, "/"))   # V x 3

  # --- Build coordinate -> voxel index hash ------------------------------
  keys      <- paste(grid[, 1L], grid[, 2L], grid[, 3L], sep = "_")
  coord_map <- new.env(hash = TRUE, parent = emptyenv(), size = V)
  for (v in seq_len(V)) coord_map[[keys[v]]] <- v

  # --- 6-connected offsets -----------------------------------------------
  offsets <- list(
    c( 1L,  0L,  0L), c(-1L,  0L,  0L),
    c( 0L,  1L,  0L), c( 0L, -1L,  0L),
    c( 0L,  0L,  1L), c( 0L,  0L, -1L)
  )

  # --- Build neighbor list and edge list ---------------------------------
  neighbor_list  <- vector("list", V)
  edge_list_rows <- vector("list", V)

  for (v in seq_len(V)) {
    gv <- grid[v, ]
    nb <- integer(0)
    for (off in offsets) {
      nb_key <- paste(gv[1L] + off[1L], gv[2L] + off[2L], gv[3L] + off[3L],
                      sep = "_")
      nb_v <- coord_map[[nb_key]]
      if (!is.null(nb_v)) nb <- c(nb, nb_v)
    }
    neighbor_list[[v]] <- nb
    new_edges <- nb[nb > v]
    if (length(new_edges) > 0L)
      edge_list_rows[[v]] <- cbind(v, new_edges)
  }

  .sira_finalise_neighbor_structures(neighbor_list, edge_list_rows, V, env)
}


# -----------------------------------------------------------------------------
# .sira_finalise_neighbor_structures
# Shared post-processing: edge list, edge index, smoothing matrix.
# Called by both neighbor builders.
# -----------------------------------------------------------------------------

#' @keywords internal
.sira_finalise_neighbor_structures <- function(neighbor_list, edge_list_rows,
                                               V, env) {
  edge_list <- do.call(rbind, edge_list_rows)
  storage.mode(edge_list) <- "integer"
  colnames(edge_list) <- c("v1", "v2")

  edge_structure <- .create_edge_index(edge_list)

  i_idx   <- edge_list[, 1L]
  j_idx   <- edge_list[, 2L]
  adj_i   <- c(i_idx, j_idx)
  adj_j   <- c(j_idx, i_idx)
  adj_mat <- Matrix::sparseMatrix(
    i = adj_i, j = adj_j, x = rep(1, length(adj_i)),
    dims = c(V, V)
  )
  degrees  <- Matrix::rowSums(adj_mat)
  deg_safe <- pmax(degrees, 1L)
  sm <- Matrix::Diagonal(x = 0.8 / deg_safe) %*% adj_mat +
    Matrix::Diagonal(x = rep(0.2, V))

  env$neighbor_list    <- neighbor_list
  env$edge_list        <- edge_list
  env$edge_structure   <- edge_structure
  env$smoothing_matrix <- sm

  if (isTRUE(env$verbose))
    message(sprintf("  Voxels: %d   |   Edges: %d   |   Avg neighbors: %.2f",
                    V, nrow(edge_list), mean(degrees)))
}


# -----------------------------------------------------------------------------
# .voxel_neighbors  (grid layout only)
# -----------------------------------------------------------------------------

#' @keywords internal
.voxel_neighbors <- function(v, d1, d2, d3) {
  i3  <- ceiling(v / (d1 * d2))
  rem <- v - (i3 - 1L) * d1 * d2
  i2  <- ceiling(rem / d1)
  i1  <- rem - (i2 - 1L) * d1

  nb <- integer(6)
  k  <- 0L
  if (i1 > 1L)  { k <- k + 1L; nb[k] <- v - 1L }
  if (i1 < d1)  { k <- k + 1L; nb[k] <- v + 1L }
  if (i2 > 1L)  { k <- k + 1L; nb[k] <- v - d1 }
  if (i2 < d2)  { k <- k + 1L; nb[k] <- v + d1 }
  if (i3 > 1L)  { k <- k + 1L; nb[k] <- v - d1 * d2 }
  if (i3 < d3)  { k <- k + 1L; nb[k] <- v + d1 * d2 }
  nb[seq_len(k)]
}


# -----------------------------------------------------------------------------
# .voxel_array_index  (grid layout only)
# Vectorised linear -> (i1, i2, i3) for use in the split operation.
# -----------------------------------------------------------------------------

#' @keywords internal
.voxel_array_index <- function(voxels, d1, d2) {
  i1 <- ifelse(voxels %% d1 != 0L, voxels %% d1, d1)
  i2 <- ceiling(((voxels - 0.1) %% (d1 * d2)) / d1)
  i3 <- ceiling(voxels / (d1 * d2))
  matrix(c(i1, i2, i3), ncol = 3L,
         dimnames = list(NULL, c("i1", "i2", "i3")))
}


# -----------------------------------------------------------------------------
# .get_voxel_coords
# Returns an |voxels| x 3 coordinate matrix for `voxels`, using either
# env$coords (coords layout) or .voxel_array_index (grid layout).
# Used by the split operation in region_ops.R.
# -----------------------------------------------------------------------------

#' @keywords internal
.get_voxel_coords <- function(voxels, env) {
  if (!is.null(env$coords)) {
    env$coords[voxels, , drop = FALSE]
  } else {
    .voxel_array_index(voxels, env$d1, env$d2)
  }
}


# -----------------------------------------------------------------------------
# .get_region_neighbors / .get_connected_components / .get_largest_...
# Unchanged — operate on neighbor_list which is the same in both layouts.
# -----------------------------------------------------------------------------

#' @keywords internal
.get_region_neighbors <- function(region_voxels, env) {
  nb <- unique(unlist(env$neighbor_list[region_voxels], use.names = FALSE))
  nb[!(nb %in% region_voxels)]
}

#' @keywords internal
.create_edge_index <- function(edge_list) {
  n_edges       <- nrow(edge_list)
  unique_voxels <- unique(c(edge_list[, 1L], edge_list[, 2L]))
  edge_index    <- new.env(hash = TRUE, parent = emptyenv(),
                           size = length(unique_voxels))
  for (i in seq_len(n_edges)) {
    v1 <- as.character(edge_list[i, 1L])
    v2 <- as.character(edge_list[i, 2L])
    if (exists(v1, envir = edge_index, inherits = FALSE)) {
      edge_index[[v1]] <- c(edge_index[[v1]], i)
    } else {
      edge_index[[v1]] <- i
    }
    if (exists(v2, envir = edge_index, inherits = FALSE)) {
      edge_index[[v2]] <- c(edge_index[[v2]], i)
    } else {
      edge_index[[v2]] <- i
    }
  }
  list(edges = edge_list, index = edge_index)
}

#' @keywords internal
.get_connected_components <- function(voxels, neighbor_list) {
  if (length(voxels) == 0L) return(list())
  voxel_to_idx <- stats::setNames(seq_along(voxels), voxels)
  visited      <- logical(length(voxels))
  components   <- list()
  for (start in seq_along(voxels)) {
    if (visited[start]) next
    queue     <- start
    component <- integer(0)
    while (length(queue) > 0L) {
      cur          <- queue[1L]
      queue        <- queue[-1L]
      if (visited[cur]) next
      visited[cur] <- TRUE
      component    <- c(component, cur)
      nb_vox <- neighbor_list[[voxels[cur]]]
      nb_idx <- voxel_to_idx[as.character(nb_vox)]
      nb_idx <- nb_idx[!is.na(nb_idx) & !visited[nb_idx]]
      queue  <- c(queue, nb_idx)
    }
    components[[length(components) + 1L]] <- voxels[component]
  }
  components
}

#' @keywords internal
.get_largest_connected_component <- function(voxels, neighbor_list) {
  if (length(voxels) <= 1L) return(voxels)
  comps <- .get_connected_components(voxels, neighbor_list)
  sizes <- lengths(comps)
  comps[[which.max(sizes)]]
}
