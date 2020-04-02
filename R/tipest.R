#Working in class on way to estimate tips. Showing simulation and how they can be useful

# library(ape)
# library(phytools)
# library(geiger)
#
#
# # -a  a
# #  b -b
# r12 <- .08
# r21 <- .06
# ntax <- 50
#
# phy <- ape::rcoal(ntax)
# q <- list(rbind(c(-r12, r12), c(r21, -r21)))
# traits <- geiger::sim.char(phy=phy, par=q, nsim=1, model="discrete", root=1)[,,1]
# taxon_to_delete <- sample(sequence(ntax), 1)
# traits_changed <- traits[-taxon_to_delete]
# phy_changed <- ape::drop.tip(phy, names(traits)[taxon_to_delete])
# model_fit <- geiger::fitDiscrete(phy=phy_changed, dat=traits_changed, model="ARD", transform="none")
#
# traits_1 <- traits
# traits_1[taxon_to_delete] <- 1
# model_fit_1 <- geiger::fitDiscrete(phy=phy, dat=traits_1, model="ARD", transform="none")
#
# traits_2 <- traits
# traits_2[taxon_to_delete] <- 2
# model_fit_2 <- geiger::fitDiscrete(phy=phy, dat=traits_2, model="ARD", transform="none")
#
# print(taxon_to_delete)
# print(rbind(traits, traits_1, traits_2))
# print(model_fit_1$opt$aicc)
# print(model_fit_2$opt$aicc)

#' Infer state at one tip
#'
#' This function estimates the best tip state (given AICc) but also gives you access to the AICc values for all the possible tip states.
#'
#' Note that this gives the best solution given the observed tips and ONLY this tip. It does not compare with other tips to reconstruct to get the best global solution.
#' @param tip_to_reconstruct Label of the tip to reconstruct
#' @param existing_tips Vector of tip states, with labels equal to tips on the tree
#' @param rate_matrix Rate matrix to pass to geiger for reconstruction
#' @param phy A phylo object (your tree)
#' @return A list with $best (the state of the best tip) and $all, a data.frame of all tip states and their AICc values
#' @export
EstimateOneTip <- function(tip_to_reconstruct, existing_tips, rate_matrix, phy) {
  tips_to_keep <- c(tip_to_reconstruct, names(existing_tips))
  phy_pruned <- ape::keep.tip(phy, tips_to_keep)
  existing_tips_with_missing <- c(existing_tips, NA)
  names(existing_tips_with_missing)[length(existing_tips_with_missing)] <- tip_to_reconstruct
  possible_states <- sort(unique(existing_tips))
  results <- data.frame(state=possible_states, AICc=NA, stringsAsFactors=FALSE)
  for (state_index in sequence(length(possible_states))) {
    filled_in_tips <- existing_tips_with_missing
    filled_in_tips[length(filled_in_tips)] <- possible_states[state_index]
    model_fit <- geiger::fitDiscrete(phy=phy, dat=filled_in_tips, model=rate_matrix, transform="none")
    results$AICc[state_index] <- model_fit$opt$aicc
  }
  best <- results$state[which.min(results$AICc)]
  return(list(best=best, all=results))
}
