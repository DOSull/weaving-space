bind_rows(
  rect11_unit$primitive,
  rect32_unit$primitive %>% sf_shift(c(750, 0), .),
  rect34_unit$primitive %>% sf_shift(c(0, -750), .),
  twill_unit$primitive %>% sf_shift(c(1050, -750), .), 
  basket_unit$primitive %>% sf_shift(c(1800, -750), .),
  this_unit$primitive %>% sf_transform(wk_affine_scale(2, 2)) %>% sf_shift(c(750, -2100), .),
) %>% plot(border = NA, main = "Some biaxial weaves")
