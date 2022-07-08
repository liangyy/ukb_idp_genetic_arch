mm <- png::readPNG('~/Documents/misc/brainxcan_figures/192915')
c1 <- c(26, 138, 43)
c2 <- c(209, 209, 209)
change_color <- function(mat, color, mask = NULL) {
  for(i in 1 : 3) {
    tmp <- mat[, ,i]
    if(is.null(mask)) tmp[tmp > 0] <- color[i] / 255
    else tmp[mask] <- color[i] / 255
    mat[, , i] <- tmp
  }
  return(mat)
}
green <- change_color(mm, c1)
gray <- change_color(mm, c2)
png::writePNG(green, '~/Documents/misc/brainxcan_figures/192915_green.png')
png::writePNG(gray, '~/Documents/misc/brainxcan_figures/192915_gray.png')

mm2 <- png::readPNG('~/Documents/misc/brainxcan_figures/273487')
add_circle <- function(mat, center, radius, color) {
  w <- ncol(mat[, , 1])
  h <- nrow(mat[, , 1])
  center <- c(round(w * center[1]), round(h * center[2]))
  radius <- w * radius
  for(k in 1 : 3) {
    tmp <- mat[, , k]
    for(i in 1 : w) {
      for(j in 1 : h) {
        dd <- sum((center[2:1] - c(j, i)) ^ 2) <= radius ^ 2
        if(dd) tmp[j, i] <- color[k] / 255
      }
    } 
    mat[, , k] <- tmp
  }
  return(mat)
}
add_line <- function(mat, from, to, width, color) {
  w <- ncol(mat[, , 1])
  h <- nrow(mat[, , 1])
  from <- c(round(w * from[1]), round(h * from[2]))
  to <- c(round(w * to[1]), round(h * to[2]))
  width <- w * width
  dx <- from[1] - to[1]
  dy <- from[2] - to[2]
  d <- sqrt(dx ^ 2 + dy ^ 2)
  slope <- dy / dx
  ddx <- width / d * dx
  ddy <- width / d * dy
  line0 <- function(x0, a, point) {
    return(point[2] - a * (point[1] - x0))
  }
  line_up <- function(x0) {
    pp <- from + c(ddx, -ddy)
    return(line0(x0, slope, pp))
  }
  line_bt <- function(x0) {
    pp <- from + c(-ddx, ddy)
    return(line0(x0, slope, pp))
  }
  for(k in 1 : 4) {
    tmp <- mat[, , k]
    for(i in from[1] : to[1]) {
      j <- c(round(line_up(i)), round(line_bt(i)))
      tmp[j[1] : j[2], i] <- ifelse(k < 4, color[k] / 255, 1)
    } 
    mat[, , k] <- tmp
  }
  return(mat)
}
blue <- c(123, 173, 255)
yellow <- c(249, 249, 98)
c3 <- c(188, 188, 188)
# br <- mm2
br <- change_color(mm2, c2, mm2[, , 1] > 0)
br2 <- br
br2 <- add_line(br2, c(0.35, 0.3), c(0.7, 0.7), 0.02, blue)
br2 <- add_line(br2, c(0.65, 0.3), c(0.3, 0.7), 0.02, blue)
br2 <- add_circle(br2, c(0.35, 0.3), 0.07, yellow)
br2 <- add_circle(br2, c(0.65, 0.3), 0.07, yellow)
br2 <- add_circle(br2, c(0.3, 0.7), 0.1, yellow)
br2 <- add_circle(br2, c(0.7, 0.7), 0.1, yellow)
png::writePNG(br2, '~/Documents/misc/brainxcan_figures/273487_br.png')
