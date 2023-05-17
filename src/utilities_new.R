egp_define_counterfactual <-
  function(data,
           som,
           cf.ids = NULL,
           fac.ids = NULL,
           compare.by = NULL,
           group.by = NULL,
           som.var = "som.unit",
           geo.vars = NULL,
           geo.kernel = "matern32",
           geo.range = NULL,
           id.var = "id",
           group.name = "group.id",
           unit.name = "cf.unit",
           assign.name = "assigned",
           assign.cat = c("counterfactual", "factual"),
           n.min = 1,
           nb.strategy = "sequential",
           deg.max = NULL,
           agg.size = NULL,
           progress = TRUE) {

  data.dt <- as.data.table(data)

  if(!id.var %in% names(data.dt)) {
    data.dt[, id.col := 1:.N, env = list(id.col = id.var)]
  }

  if(is.list(group.by)) {
    group.by.c <- unique(do.call(c, group.by))
  }
  if(!is.list(group.by)) {
    group.by.c <- group.by
    group.by <- list(group.by)
  }

  vars.sel <- c(id.var, compare.by, group.by.c, som.var, geo.vars)

  data.dt <- data.dt[, ..vars.sel]

  if(is.null(fac.ids) & is.null(cf.ids)) {
    stop("Either `fac.ids` or `cf.ids` must be provided")
  }
  if(!is.null(fac.ids) & is.null(cf.ids)) {
    cf.ids <-
      data.dt[!id.col %in% fac.ids,
              id.col,
              env = list(id.col = id.var)]
  }
  if(is.null(fac.ids) & !is.null(cf.ids)) {
    fac.ids <-
      data.dt[!id.col %in% cf.ids,
              id.col,
              env = list(id.col = id.var)]
  }

  # Reference SOM units

  data.dt <-
    data.dt[.(id.col = unique(c(fac.ids, cf.ids))),
            on = id.var,
            nomatch = NULL,
            env = list(id.col = id.var)]
  setkeyv(data.dt, id.var)

  bmu.ref <- data.dt[.(cf.ids), on = id.var]

  if(is.null(compare.by)) {

    n.bmu <- integer(nrow(som$grid.nb))
    n.ref <- bmu.ref[order(som.col),
                     .(n = .N),
                     by = som.col,
                     env = list(som.col = som.var)]
    n.bmu[n.ref[, som.col, env = list(som.col = som.var)]] <- n.ref$n
    n.bmu[is.na(n.bmu)] <- 0


    if(nb.strategy == "sequential") {

      dist = som$unit.dist
      n = n.bmu
      n.min = n.min
      deg.max = deg.max

      nbh <-
        .nb_sequential(dist = som$unit.dist,
                       n = n.bmu,
                       n.min = n.min,
                       deg.max = deg.max)
    }
    if(nb.strategy == "expand") {
      nbh <-
        .nb_expand(A = som$grid.nb,
                   n = n.bmu,
                   n.min = n.min,
                   deg.max = deg.max)
    }

    cf.bmu.dt <- 
      data.dt[.(fac.ids),
              on = id.var
              ][,
                .(id.col,
                  cf.col = nbh$neighbourhood[som.col]),
                env = list(id.col = id.var,
                           cf.col = unit.name,
                           som.col = som.var)]
    cf.bmu.dt[,
              .n := as.integer(lapply(cf.col, length)),
              env = list(cf.col = unit.name)]

    cf.ids.dt <-
      data.dt[.(cf.ids),
              on = id.var
              ][order(som.col, id.col),
                .(id.col = list(id.col)),
                by = (cf.col = som.col),
                env = list(id.col = id.var,
                           som.col = som.var,
                           cf.col = unit.name)]

  } else {


    bmus <- 1:nrow(som$grid.nb)

    bmu.by <-
      data.dt[,
              .(som.col = 1:nrow(som$grid.nb)),
              by = compare.by,
              env = list(som.col = som.var)]

    comp <-
      data.dt[, .(.n = .N), by = c(compare.by)]
    comp[, .cfgrp := 1:.N] 

    comp.bmu.cf <-
      data.dt[.(cf.ids), on = id.var
              ][, .(.n = .N), by = c(som.var, compare.by)] |>
      merge(bmu.by, by = c(som.var, compare.by), all = TRUE)
 
    nbh.l <- list()
    cf.ids.l <- list()
    for(i in 1:nrow(comp)) {
      
      n.ref <- comp.bmu.cf[comp[i, compare.by, with = FALSE],
                           on = compare.by
                           ][order(som.col),
                             .(som.col, .n),
                             env = list(som.col = som.var)]

      n.bmu <- integer(nrow(som$grid.nb))
      n.bmu[n.ref[, som.col, env = list(som.col = som.var)]] <- n.ref$.n
      n.bmu[is.na(n.bmu)] <- 0

      if(nb.strategy == "expand") {
        nbh.l[[i]] <-
          .nb_expand(A = som$grid.nb,
                     n = n.bmu,
                     n.min = n.min,
                     deg.max = deg.max)$neighbourhood
      }
      if(nb.strategy == "expand") {
        nbh.l[[i]] <-
          .nb_sequential(dist = som$unit.dist,
                         n = n.bmu,
                         n.min = n.min,
                         deg.max = deg.max)$neighbourhood
      }
      
      cf.ids.l[[i]] <-
        data.dt[.(cf.ids),
                on = id.var
                ][comp[i, compare.by, with = FALSE],
                  on = compare.by
                  ][order(som.col, id.col),
                    .(id.col = list(id.col)),
                    by = .(cf.col = som.col),
                    env = list(id.col = id.var,
                               som.col = som.var,
                               cf.col = unit.name)]

    }

    comp[, `:=`(.nbh = nbh.l)]

    comp.cf <-
      comp[,
           .(cf.col = unlist(.nbh, recursive = FALSE),
             som.col = bmus),
           by = c(".cfgrp", compare.by),
           env = list(som.col = som.var,
                      cf.col = unit.name)
           ]

    cf.bmu.dt <-
      cbind(
            data.dt[.(fac.ids),
                    on = id.var
                    ][order(id.col),
                      env = list(id.col = id.var)
                      ][,
                        c(id.var, compare.by), with = FALSE],
            comp.cf[data.dt[.(fac.ids),
                            on = id.var
                            ][order(id.col),
                              env = list(id.col = id.var)
                              ][,
                                c(compare.by, som.var),
                                with = FALSE],
                                .(cf.col),
                                on = c(compare.by, som.var),
                                env = list(cf.col = unit.name)])
 
    cf.ids.dt <- 
      rbindlist(cf.ids.l, idcol = ".cfgrp") |>
      merge(comp[, -c(".n", ".nbh")], by = ".cfgrp") |>
      setcolorder(c(compare.by, unit.name, id.var)) |>
      DT(, -".cfgrp")

  }

  cf.ids.dt[,
            .n := unlist(lapply(id.col, length)),
            env = list(id.col = id.var)]


  data.cf <- data.dt[.(cf.ids), on = id.var ]
  data.cf[,
          assign.col := assign.cat[1],
          env = list(assign.col = assign.name)]

  data.fac <- data.dt[.(fac.ids), on = id.var]
  data.fac[,
           assign.col := assign.cat[2],
           env = list(assign.col = assign.name)]
  
  data.ret <- rbind(data.fac, data.cf)
  setorderv(data.ret, id.var)

  # Groups

  id.names <- paste(id.var, assign.cat, sep = ".")
  geo.names.cf <- paste(geo.vars, assign.cat[1], sep = ".")
  geo.names.fac <- paste(geo.vars, assign.cat[2], sep = ".")

  groups.l <- list()
  for(i in seq_along(group.by)){
    groups.l[[i]] <-
      .ids_by_group(data.dt[.(fac.ids)],
                    id.var = id.var,
                    group.vars = group.by[[i]],
                    add.label = FALSE)
  }
  groups <- rbindlist(groups.l, fill = TRUE)
  setnames(groups, id.var, id.names[2])
  groups[,
         group.col := 1:nrow(groups),
         env = list(group.col = group.name)]
  setorderv(groups, group.name)

  groups.long <-
    groups[,
           .(id.fac = as.integer(unlist(cell.factual))),
           by = group.name,
           env = list(id.fac = id.names[2])]
  groups.long[, .n.fac := .N, by = group.name]
  setkeyv(groups.long, id.names[2])

  # groups.n <-
  #   groups.long[, .(.n.fac = .N), by = group.name]
  # setkey(groups.n, .n.fac)


  units.fac <-
    cf.bmu.dt[,
               .(cf.col = as.integer(unlist(cf.col))),
               by = id.var,
               env = list(cf.col = unit.name)
               ][,
                 .(id.fac = id.col,
                   cf.col),
                 env = list(id.fac = id.names[2],
                            id.col = id.var,
                            cf.col = unit.name)] |>
    merge(SJ(cf.ids.dt[,
                       .(cf.col, .n),
                       env = list(cf.col = unit.name)]),
          by = unit.name,
          sort = FALSE,
          all.x = TRUE,
          all.y = FALSE) |>
    merge(data.fac[,
                   .(id.fac = id.col,
                     fac.x = geo.x,
                     fac.y = geo.y),
                   env = list(id.fac = id.names[2],
                              id.col = id.var,
                              fac.x = geo.names.fac[1],
                              geo.x = geo.vars[1],
                              fac.y = geo.names.fac[2],
                              geo.y = geo.vars[2])],
          by = id.names[2]) |>
    na.omit(".n")


  # REMOVE later
  if(!all(data.fac$cell %in% units.fac$cell.factual)) warnings("not all obs assigned")

  units.cf <-
    cf.ids.dt[,
               .(id.cf = as.integer(unlist(id.col))),
               by = cf.col,
               env = list(id.cf = id.names[1],
                          id.col = id.var,
                          cf.col = unit.name)
               ] |>
    merge(data.cf[,
                  .(id.cf = id.col,
                    cf.x = geo.x,
                    cf.y = geo.y),
                  env = list(id.cf = id.names[1],
                             id.col = id.var,
                             cf.x = geo.names.cf[1],
                             geo.x = geo.vars[1],
                             cf.y = geo.names.cf[2],
                             geo.y = geo.vars[2])],
          by = id.names[1])


  if(is.null(agg.size)) agg.size <- sum(units.fac$.n)

  units.fac[order(-.n), .nc := cumsum(.n)]
  units.fac[, .agg.id := ceiling(.nc / agg.size)]
  units.fac <- units.fac[, -c(".n", ".nc")]
  agg.ids <- na.omit(unique(units.fac[order(.agg.id), .agg.id]))

  if(!is.null(geo.vars)) {
    range.sam <- 1e6
    if(is.null(geo.range)) {
      if(nrow(data.dt) > range.sam) {
        pts <- st_multipoint(x = as.matrix(data.dt[sample(1:.N, range.sam, replace = TRUE),
                                                   ..geo.vars]),
                             dim = "XY")
      } else {
        pts <-
          st_multipoint(x = as.matrix(data.dt[, ..geo.vars]),
                             dim = "XY")
      }
      pts.bb <- st_bbox(st_minimum_bounding_circle(pts))
      geo.range <- pts.bb[["xmax"]] - pts.bb[["xmin"]]
    }
  } else {
    geo.kernel = "nodist"
    geo.range <- 1
  }

  dist.fun <- function(x) sqrt(rowSums((x[, 1:2] - x[, 3:4])^2))

  geo.fun <-
    switch(geo.kernel,
           "matern12" = function(x) {exp(-x)},
           "matern32" = function(x) {
                          sqrt3 <- sqrt(3)
                          (1 + sqrt3 * x) * exp(-sqrt3 * x)},
           "matern52" = function(x) {
                          sqrt5 <- sqrt(5)
                          (1 + sqrt5 * x + sqrt5/3 * x^2) * exp(-sqrt5 * x)},
          )


  if(length(agg.ids) <= 1) progress <- FALSE

  if(progress) {
    prog <- txtProgressBar(min = 0,
                           max = length(agg.ids),
                           char = "=", width = NA, title = "Progress", style = 3)
  }

  weights.cf.l <- list()

  for(i in seq_along(agg.ids)) {

    if(!is.null(geo.vars)) {

      dist.obs <-
        units.cf[SJ(units.fac[J(agg.ids[i]),
                    on = ".agg.id"]),
                 on = unit.name,
                 allow.cartesian = TRUE
                 ][,
                   .(fac.id,
                     cf.id,
                     .dist = dist.fun(as.matrix(.SD)) / geo.range),
                   .SDcols = c(geo.names.fac, geo.names.cf),
                   env = list(fac.id = id.names[2],
                              cf.id = id.names[1])]
    dist.obs[,
             .geosim := geo.fun(.dist)]

    } else {

      dist.obs <-
        units.cf[SJ(units.fac[J(agg.ids[i]),
                    on = ".agg.id"]),
                 on = unit.name,
                 allow.cartesian = TRUE
                 ][,
                   .(fac.id,
                     cf.id,
                     .geosim = 1),
                   env = list(fac.id = id.names[2],
                              cf.id = id.names[1])]


    }

    dist.obs[,
             .w := .geosim / sum(.geosim),
             by = cell.factual]
    setkeyv(dist.obs, id.names[2])

    weights <- groups.long[dist.obs, on = id.names[2], allow.cartesian = TRUE]

    weights.cf.l[[i]] <-
      weights[,
              .(.w = sum(.w)),
              by = c(group.name, id.names[1])]

    if(progress) {
      setTxtProgressBar(prog, i)
    }

  }

  if(progress) {
    close(prog)
  }

  weights.cf <-
    rbindlist(weights.cf.l) |>
    DT(, .(.w = sum(.w)),
       by = c(group.name, id.names[1]))

  weights.cf[, .w := .w/sum(.w), by = group.name]

  groups.cf <-
    weights.cf[order(cf.id),
               .(cf.id = list(cf.id),
                 .w = list(.w)),
               by = group.name,
               env = list(cf.id = id.names[1])
               ][order(group.col),
                 env = list(group.col = group.name)]

  groups.comb <- merge(groups, groups.cf)
  setcolorder(groups.comb, c(group.name, group.by.c, id.names[2], id.names[1], ".w")) 

  counterfactual <- list(data = data.ret,
                         groups = groups.comb,
                         assignment = cf.bmu.dt[, -".n"],
                         units = cf.ids.dt[, -".n"],
                         compare.by = compare.by,
                         group.by = group.by,
                         group.by.c = group.by.c,
                         id.var = id.var,
                         som.var = som.var,
                         geo.vars = geo.vars,
                         geo.kernel = geo.kernel,
                         geo.range = geo.range,
                         group.var = group.name,
                         unit.var = unit.name,
                         assign.var = assign.name,
                         assign.cat = assign.cat)

  class(counterfactual) <- c(class(counterfactual), "egp_cf_def")

  return(counterfactual)
}

egp_evaluate_factual <- function(predictions,
                                 cf.def,
                                 name = "factual",
                                 group.eval = NULL,
                                 pred.var = NULL,
                                 draw.chunk = NULL,
                                 agg.size = NULL,
                                 parallel = NULL,
                                 progress = TRUE,
                                 ...
                                 ){

  predictions.dt <- as.data.table(predictions)

  id.var <- cf.def$id.var
  id.fac <- paste(id.var, cf.def$assign.cat[2], sep = ".")
  id.cf <- paste(id.var, cf.def$assign.cat[1], sep = ".")
  group.var <- cf.def$group.var

  if(is.null(pred.var)) {
    pred.sel <- which(!names(predictions.dt) %in% c(id.var, ".draw"))[1]
    pred.var <- names(predictions.dt)[pred.sel]
  }

  if(is.null(group.eval)) {
    group.eval <- cf.def$groups[[group.var]]
  }

  factual <-
    .aggregate_variables(predictions.dt,
               agg.fun = mean,
               ids = cf.def$groups[.(group.eval),
                                   id.col,
                                   on = group.var,
                                   env = list(id.col = id.fac)],
               pred.var = pred.var,
               draw.var = ".draw",
               id.var = id.var,
               agg.name = name,
               group.name = group.var,
               draw.chunk = draw.chunk,
               agg.size = agg.size,
               parallel = parallel,
               progress = progress)

  factual[,
          group.col := group.eval[group.col],
          env = list(group.col = group.var)]
  factual.aug <-
    merge(cf.def$groups[, -c(id.fac, id.cf, ".w"), with = FALSE], factual, all.x = FALSE)

  setkeyv(factual.aug, group.var)
  setindexv(factual.aug, ".draw")
  return(factual.aug)
}


egp_evaluate_counterfactual <- function(predictions,
                                        cf.def,
                                        name = "counterfactual",
                                        group.eval = NULL,
                                        pred.var = NULL,
                                        draw.chunk = NULL,
                                        agg.size = NULL,
                                        parallel = NULL,
                                        progress = TRUE,
                                        ...
                                        ){

  predictions.dt <- as.data.table(predictions)

  id.var <- cf.def$id.var
  id.fac <- paste(id.var, cf.def$assign.cat[2], sep = ".")
  id.cf <- paste(id.var, cf.def$assign.cat[1], sep = ".")
  group.var <- cf.def$group.var

  if(is.null(pred.var)) {
    pred.sel <- which(!names(predictions.dt) %in% c(id.var, ".draw"))[1]
    pred.var <- names(predictions.dt)[pred.sel]
  }

  if(is.null(group.eval)) {
    group.eval <- cf.def$groups[[group.var]]
  }

  counterfactual <-
    .aggregate_variables(predictions.dt,
               agg.fun = sum,
               ids = cf.def$groups[.(group.eval),
                                   id.col,
                                   on = group.var,
                                   env = list(id.col = id.cf)],
               weights = cf.def$groups[.(group.eval),
                                       .w,
                                       on = group.var],
               pred.var = pred.var,
               draw.var = ".draw",
               id.var = id.var,
               agg.name = name,
               group.name = group.var,
               draw.chunk = draw.chunk,
               agg.size = agg.size,
               parallel = parallel,
               progress = progress)

  counterfactual[,
                 group.col := group.eval[group.col],
                 env = list(group.col = group.var)]
  counterfactual.aug <-
    merge(cf.def$groups[, -c(id.fac, id.cf, ".w"), with = FALSE], counterfactual, all.x = FALSE)

  setkeyv(counterfactual.aug, group.var)
  setindexv(counterfactual.aug, ".draw")
  return(counterfactual.aug)
}


.aggregate_variables.data.table <- 
  function(predictions,
           agg.fun = mean,
           trans.fun = NULL,
           ids = NULL,
           weights = NULL,
           pred.var = "predicted",
           draw.var = ".draw",
           id.var = "id",
           agg.name = "aggregated",
           group.name = "group.id",
           draw.ids = NULL,
           draw.chunk = NULL,
           agg.size = NULL,
           parallel = NULL,
           progress = TRUE,
           ...) {
  predictions.dt <- as.data.table(predictions)
  setkeyv(predictions.dt, id.var)
  setindexv(predictions.dt, draw.var)
  if(is.numeric(parallel)) {
    dt.threads.old <- getDTthreads()
    arrow.threads.old <- cpu_count()
    setDTthreads(parallel)
    set_cpu_count(parallel)
  }
  if(is.null(draw.ids)) {
    draw.ids <- predictions.dt[, unique(draw.col), env = list(draw.col = draw.var)]
  }
  if(is.null(draw.chunk)) {
    draw.chunk <- length(draw.ids)
  }
  if(is.null(ids)) {
    ids <- list(predictions.dt[, unique(id.col), env = list(id.col = id.var)])
  }
  if(is.null(weights)) {
    weights <- lapply(ids, \(x) rep(1, length(x)))
  }
  draw.chunks <- chunk_seq(1, length(draw.ids), draw.chunk)
  if(is.null(agg.size))
    agg.size <- min(unlist(lapply(ids, length)))
  ids.dt <-
    list(1:length(ids), ids, weights) |>
    as.data.table() |>
    setnames(c(group.name, id.var, ".w"))
  setorderv(ids.dt, group.name)
  ids.dt[,N := as.numeric(unname(unlist(lapply(ids, length))))]
  ids.dt[order(-N), Nc := cumsum(N)]
  ids.dt[,
         `:=`(aggregate = ifelse(N < agg.size, TRUE, FALSE))
         ][aggregate == TRUE, 
           agg.id := ceiling(Nc / agg.size)]
  setkeyv(ids.dt, group.name)
  agg.ids <- na.omit(unique(ids.dt[order(agg.id), agg.id]))
  single.ids <- ids.dt[aggregate == FALSE,
                       group.col,
                       env = list(group.col = group.name)]
  if(all(ids.dt$N == 0)) return(NULL)
  ids.dt.l <-
    ids.dt[order(-N),
           .(id.col = unlist(id.col),
             .w = unlist(.w)),
           by = c(group.name, "agg.id"),
           env = list(id.col = id.var)]
  setkeyv(ids.dt.l, id.var)
  setindexv(ids.dt.l, list(group.name, "agg.id"))
  if(!is.null(trans.fun)) {
    predictions.dt[,
                   pred.col := trans(pred.col),
                   env = list(pred.col = pred.var,
                              trans = trans.fun)]
  }
  if(progress) {
    prog <- txtProgressBar(min = 0, max = max(ids.dt$Nc) * length(draw.chunks$from), initial = 0,
                           char = "=", width = NA, title = "Progress", style = 3)
    prog.counter <- 0
  }
  draws.agg.l <- list()
  for(i in seq_along(draw.chunks$from)) {
    draws.sum <- draw.ids[draw.chunks$from[i]:draw.chunks$to[i]]
    predictions.draws <- predictions.dt[.(draws.sum),
                                        on = draw.var,
                                        nomatch = NULL]
    setkeyv(predictions.draws, id.var)
    single.l <- list()
    for(j in seq_along(single.ids)) {
      ids.sum <- ids.dt.l[.(single.ids[j]), on = group.name]
      single.l[[j]] <-
        predictions.draws[.(ids.sum),
                          nomatch = NULL,
                          on = id.var
                          ][order(draw.col, group.col),
                            .(agg.col = agg.fun(.w*pred.col)),
                            by = c(draw.var, group.name),
                            env = list(group.col = group.name,
                                       draw.col = draw.var,
                                       pred.col = pred.var,
                                       agg.col = agg.name,
                                       agg = agg.fun)]
      if(progress) {
        prog.counter <- prog.counter + ids.dt[group.col == single.ids[j],
                                              N,
                                              env = list(group.col = group.name)]
        setTxtProgressBar(prog, prog.counter)
      }
    }
    agg.l <- list()
    for(k in seq_along(agg.ids)) {
      match.ids.groups <- ids.dt.l[.(agg.ids[k]), on = "agg.id"]
      agg.l[[k]] <-
        predictions.draws[match.ids.groups,
                          nomatch = NULL,
                          on = id.var,
                          allow.cartesian = TRUE
                          ][order(draw.col, group.col),
                            .(agg.col = agg.fun(.w*pred.col)),
                            by = c(draw.var, group.name),
                            env = list(group.col = group.name,
                                       draw.col = draw.var,
                                       pred.col = pred.var,
                                       agg.col = agg.name,
                                       agg = agg.fun)]
      if(progress) {
        prog.counter <- prog.counter + sum(ids.dt[agg.id == agg.ids[k], N])
        setTxtProgressBar(prog, prog.counter)
      }
    }
    draws.agg.l[[i]] <- rbindlist(c(single.l, agg.l), fill = TRUE)
  }
  draws.agg <- rbindlist(draws.agg.l, fill = TRUE)
  setorderv(draws.agg, c(draw.var, group.name))
  if(progress) close(prog)
  if(is.numeric(parallel)) {
    setDTthreads(dt.threads.old)
  }
  return(draws.agg)
}

