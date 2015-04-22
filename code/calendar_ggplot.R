

gg.calendar <- function(df) {
  require(ggplot2)
  require(lubridate)
  wom <- function(date) { # week-of-month
    first <- wday(as.Date(paste(year(date),month(date),1,sep="-")))
    return((mday(date)+(first-2)) %/% 7+1)
  }
  df$month <- month(df$dates)
  df$day   <- mday(df$dates)
  
  rng   <- range(df$dates)
  rng   <- as.Date(paste(year(rng),month(rng),1,sep="-"))
  start <- rng[1]
  end   <- rng[2]
  month(end) <- month(end)+1
  day(end)   <- day(end)  -1
  
  cal <- data.frame(dates=seq(start,end,by="day"))
  cal$year  <- year(cal$date)
  cal$month <- month(cal$date)
  cal$cmonth<- month(cal$date,label=T)
  cal$day   <- mday(cal$date)
  cal$cdow  <- wday(cal$date,label=T)
  cal$dow   <- wday(cal$date)
  cal$week  <- wom(cal$date)
  
  cal        <- merge(cal,df[,c("dates","counts")],all.x=T)
  
  ggplot(cal, aes(x=cdow,y=-week))+
    geom_tile(aes(fill=counts,colour="grey50"))+
    geom_text(aes(label=day),size=3,colour="grey20")+
    facet_wrap(~cmonth, ncol=3)+
    scale_fill_gradient2(low = "white", mid="blu", high = "red", na.value="white") + #change color here for heat map
    scale_color_manual(guide=F,values="grey50")+
    scale_x_discrete(labels=c("S","M","T","W","Th","F","S"))+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(panel.grid=element_blank())+
    labs(x="",y="")+
    coord_fixed()
  
#   ggplot(cal, aes(x=cdow,y=-week))+
#     geom_tile(aes(fill=cmonth,colour="grey50"))+
#     geom_text(aes(label=day),size=3,colour="grey20")+
#     facet_wrap(~cmonth, ncol=3)+
#     scale_fill_discrete()+
#     scale_color_manual(guide=F,values="grey50")+
#     scale_x_discrete(labels=c("S","M","T","W","Th","F","S"))+
#     theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
#     theme(panel.grid=element_blank())+
#     labs(x="",y="")+
#     coord_fixed()
  
}

generateData4cal <- function()
{
  set.seed(42)
  dates <- seq(as.Date("2012/01/01"), as.Date("2012/6/30"), by = "1 day")
  counts <- 1:length(dates)
  filterField <- sample(1:42,length(dates),replace=T)
  df <- data.frame(dates, counts, filterField)
  
  return(df)
}