library(ape)

coalescent_rate <- function(k, pop_size) {
    0.25 * k * (k - 1) / pop_size
}

leaf_node <- function(name, time) {
    list(type = "leaf", name = name, time = time)
}

branch_node <- function(name, children, time) {
    list(type = "branch",
         name = name,
         children = children,
         time = time,
         lengths = c(time - children[[1]]$time,
                     time - children[[2]]$time))
}

binary_parent <- function(child1, child2, time) {
    parent_name <- paste(child1$name, child2$name, sep = "-")
    branch_node(parent_name, list(child1, child2), time)
}

current_time <- 0
sampled_population <- list(leaf_node("beth", current_time),
                           leaf_node("gerry", current_time),
                           leaf_node("morty", current_time),
                           leaf_node("summer", current_time))
population_size <- 100

k <- +Inf

while (k > 2) {
    k <- length(sampled_population)
    coalescent_time <- rexp(1, coalescent_rate(k, population_size))
    current_time <- current_time + coalescent_time
    ixs <- sample.int(k, size = 2)
    parent_node <- binary_parent(sampled_population[[ixs[1]]], sampled_population[[ixs[2]]], current_time)
    if (k > 2) {
        sampled_population <- c(list(parent_node), sampled_population[-ixs])
    } else {
        sampled_population <- parent_node
    }
}

newick_helper <- function(node) {
    if (node$type == "leaf") {
        node$name
    } else if (node$type == "branch") {
        sprintf("(%s:%f,%s:%f)%s",
                newick_helper(node$children[[1]]),
                node$lengths[1],
                newick_helper(node$children[[2]]),
                node$lengths[2],
                node$name)
    }
}

newick <- function(node) {
    sprintf("%s;", newick_helper(node))
}

demo_tree <- read.tree(text = newick(sampled_population))

plot(demo_tree)
