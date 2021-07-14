dt_quantile_ramp <- function(df, cols){
    dt <- datatable(df)
    colRamp <- colorRamp(c("lightblue","white","red"))
    for(column in cols){
        x <- na.omit(df[[column]])
        brks <- quantile(x, probs = seq(.05, .95, .005))
        RGB <- colRamp(c(0, (brks-min(x))/(max(x)-min(x))))
        clrs <- apply(RGB, 1, function(rgb){
          sprintf("rgb(%s)", toString(round(rgb,0)))
        })
        dt <- dt %>% 
          formatStyle(column, backgroundColor = styleInterval(brks, clrs))
    }
    return(dt)
}
