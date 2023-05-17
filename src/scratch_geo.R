source("utilities.R")
source("utilities_new.R")

sam.frac = 0.01
ls1 <- readRDS("../landscapes/imbalance_high/data/0001.rds")
par <- readRDS("../landscapes/imbalance_high/parameters.rds")
ls <- ls1$landscape
set.seed(par$seed[1]) 
sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
ls <- na.omit(ls[sam,])


som <- egp_som(data, x.dim = 25, vars = c("z1", "z2", "z3", "z4"))

data[,
     c("som.unit", "som.x", "som.y") :=
       get_bmu(som, coord = TRUE, list = TRUE)]


data[, type.o := factor(as.character(type), levels = c("reference", "treatment"), ordered = TRUE)]

egp.k.geo <- 250
egp.max.knots.geo <- 2500
egp.k.som <- 250
egp.max.knots.som <- 250
n.threads <- c(2,1)
mod.egp <- bam(response ~
               s(x, y, bs = "gp", k = egp.k.geo,
                 xt = list(max.knots = egp.max.knots.geo)) +
               s(x, y, by = type.o, bs = "gp", k = egp.k.geo,
                 xt = list(max.knots = egp.max.knots.geo)) +
               s(som.x, som.y, bs = "gp", k = egp.k.som,
                 xt = list(max.knots = egp.max.knots.som)),
               data = data,
               select = TRUE,
               discrete = TRUE,
               nthreads = n.threads
)

mid <- data[type == "treatment" & som.unit == 15][1, cell]
cdat <- data[cell == mid | type == "reference", ]

egp.pred <-
  egp_posterior_predict(
                        model = mod.egp,
                        data = data,
                        id.var = "cell",
                        posterior = egp.post)


source("utilities_new.R")

egp.def <-
  egp_define_counterfactual(
                            data = data,
                            som = som,
                            fac.ids = data[type == "treatment", cell],
                            cf.ids = data[type == "reference", cell],
                            nb.strategy = "sequential",
                            group.by = list(NULL, "poly"),
                            som.var = "som.unit",
                            geo.vars = mod.geo.vars,
                            id.var = "cell",
                            deg.max = NULL,
                            agg.size = 1e6,
                            n.min = 1)

egp.post <- egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")
egp.fac <- egp_evaluate_factual(egp.pred, egp.def)
egp.count <- egp_evaluate_counterfactual(egp.pred, egp.def)
egp.mar <- egp_marginal(egp.fac, egp.count)
egp.mar[group.id == 1, mean(marginal)]

str(as.matrix(cdat[cell == mid, ..geo.vars])[, 1:2, drop = FALSE])

res.egp <- readRDS("../results/imbalance_high/egp_som25.rds")$marginal
res.egp[landscape == 1]


fac.ids = data[type == "treatment", cell]
cf.ids = data[type == "reference", cell]
compare.by = NULL
group.by = list(NULL, "poly")
som.var = "som.unit"
id.var = "cell"
group.name = "group.id"
unit.name = "cf.unit"
assign.name = "assigned"
assign.cat = c("counterfactual", "factual")
n.min = 1
nb.strategy = "sequential"
deg.max = NULL
geo.vars = c("x", "y")
# geo.vars = NULL
geo.kernel = "matern32"
geo.range = NULL
agg.size = 1e3

egp_define_counterfactual2 <-
  function(data,
           som,
           cf.ids = NULL,
           fac.ids = NULL,
           compare.by = NULL,
           group.by = NULL,
           som.var = "som.unit",
           geo.vars = c("x", "y"),
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
           progress = TRUE) {

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

  data.dt <- data.dt[.(id.col = unique(c(fac.ids, cf.ids))),
                     on = id.var,
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
          by = id.names[2])



  # REMOVE later
  if(!all(data.fac$cell %in% units.fac$cell.factual)) stop("not all obs assigned")

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

  # PROGRESS BAR here

  system.time({

  weights.cf.l <- list()

  for(i in seq_along(agg.ids)) {

    if(!is.null(geo.vars)) {

      dist.obs <-
        units.cf[SJ(units.fac[J(agg.ids[i]),
                    on = ".agg.id"]),
                 on = unit.name
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
                 on = unit.name
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

    weights <- groups.long[dist.obs, on = id.names[2]]

    weights.cf.l[[i]] <-
      weights[,
              .(.w = sum(.w)),
              by = c(group.name, id.names[1])]

  }

  }) 

  weights.cf <-
    rbindlist(weights.cf.l) |>
    DT(, .(.w = sum(.w)), by = c(group.name, id.names[1]))
  weights.cf[, .w := .w/sum(.w), by = group.name]

  groups.cf <-
    weights.cf[,
               .(cf.id = list(cf.id),
                 .w = list(.w)),
               by = group.name,
               env = list(cf.id = id.names[1])
               ][order(group.col),
                 env = list(group.col = group.name)]

  groups.comb <-
    cbind(groups,
          groups.cf[, -c(group.name), with = FALSE])
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

