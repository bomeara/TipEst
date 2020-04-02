test_that("return of best is a single state", {
  r12 <- .08
  r21 <- .06
  ntax <- 50
  q <- list(rbind(c(-r12, r12), c(r21, -r21)))
  phy <- ape::rcoal(ntax)

  traits <- geiger::sim.char(phy=phy, par=q, nsim=1, model="discrete", root=1)[,,1]
  while(length(unique(traits))<2) {
    traits <- geiger::sim.char(phy=phy, par=q, nsim=1, model="discrete", root=1)[,,1]
  }
  taxon_number_to_delete <- sample(sequence(ntax), 1)
  traits_changed <- traits[-taxon_number_to_delete]
  taxon_name_to_recover <-names(traits)[taxon_number_to_delete]

  result <- EstimateOneTip(tip_to_reconstruct=taxon_name_to_recover, existing_tips=traits_changed, rate_matrix="ARD", phy=phy)


  expect_equal(class(result$best), "numeric")
  expect_equal(nrow(result$all), 2)
})
