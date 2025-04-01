### Стиль ВШЭ для построения графиков в ggplot2 


# Libraries ---------------------------------------------------------------


library(showtext) # for fonts 
library(ggplot2)



# Colours  -----------------------------------------------------------------


# палитра цветов - из шаблона презентации ВШЭ
hse_colours <- c('#102D69', '#1656A6', '#BFBFBF', '#E61F3D', '#11A0D7', '#FFD746', '#029C63', '#EB681F')

# Fonts -------------------------------------------------------------------

#extrafont::font_import() # нужно запускать один раз. Либо при добавлении новых шрифтов 
#extrafont::loadfonts(device="win")

# Примечание: пакет extrafont работает только со шрифтам .ttf. Для .otf - нужно использовать пакет showtext. 

sysfonts::font_add(
  family = "HSE",
  regular = "../Technical/fonts/HSESans-Regular.otf",
  bold = "../Technical/fonts/HSESans-Bold.otf",
  italic = "../Technical/fonts/HSESans-Italic.otf",
  bolditalic = "../Technical/fonts/HSESans-SemiBold.otf"
)
sysfonts::font_add(
  family = "Times New Roman",
  regular = "../Technical/fonts/Times-New-Roman.otf"
)
font_families() 
hse_font <- 'HSE'
showtext_auto()

# стиль
theme_hse <- function(grid_lines = "vertical", legend.title = 'blank'){
  
  color.background = "#FFFFFF"
  color.grid.major = "#d9d9d9"
  color.text = "#181818"
  color.axis = "#181818"
  color.background.legend = "#FFFFFF"
  
  if(grid_lines == "vertical"){
    grid.major.x = element_blank()
    grid.major.y = element_line()
  } else {
    grid.major.x = element_line()
    grid.major.y = element_blank()
  }
  
  theme_hse <- theme_bw(base_size=11, base_family = 'Times New Roman') +
    theme(
      
      panel.background=element_rect(fill=color.background, color=color.background),
      plot.background=element_rect(fill=color.background, color=color.background),
      
      panel.border=element_rect(color=color.background),
      panel.grid.major=element_line(color=color.grid.major,size=.20),
      panel.grid.major.x=grid.major.x,
      panel.grid.major.y=grid.major.y,
      panel.grid.minor=element_blank(),
      
      axis.ticks.y=element_line(color=color.axis),
      axis.ticks.x=element_line(color=color.axis),
      axis.line.x = element_line(color=color.axis),
      axis.line.y = element_line(color=color.axis),
      
      legend.background = element_rect(fill=color.background.legend),
      legend.key = element_rect(fill=color.background, color=color.background),
      legend.text = element_text(size=rel(.9),color=color.text),
      legend.position = c(0.6, 0.95),
      legend.direction="vertical",
      
      plot.title=element_text(color=color.text, size=rel(1.3), face = "bold"),
      
      axis.text.x=element_text(size=rel(1.0),color=color.text),
      axis.text.y=element_text(size=rel(1.0),color=color.text),
      axis.title.x=element_text(size=rel(1.0),color=color.text, vjust=0),
      axis.title.y=element_text(size=rel(1.0),color=color.text, 
                                angle = 0
                                # ,margin = margin(r = -4.5, l = 0, unit = 'cm') # параметр margin r меняет положение надписи по оси Y
                                ), 
      
      plot.margin = margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm")
    )
  if (legend.title == 'blank') {
    theme_hse <- theme_hse + theme(legend.title = element_blank())
  }
  theme_hse
}


# Пример графика  -----------------------------------------------------------------


# library(xts)
# library(pdfetch) #FRED data
# 
# data_pc <- pdfetch_FRED(c("UNRATE", "CPIAUCSL"))
# # to quarterly frequency
# data_pc <- to.period(data_pc, period = "quarter", OHLC = FALSE)
# #View(data_pc)
# 
# #Transformations
# data_pc$l_cpi <- log(data_pc$CPIAUCSL)
# data_pc$unrate <- (data_pc$UNRATE)
# data_pc$inflation <- 100*diff(data_pc$l_cpi, 4)
# 
# #Add recession bars
# recessions.df = read.table(textConnection(
#   "Peak, Trough
#   1957-08-01, 1958-04-01
#   1960-04-01, 1961-02-01
#   1969-12-01, 1970-11-01
#   1973-11-01, 1975-03-01
#   1980-01-01, 1980-07-01
#   1981-07-01, 1982-11-01
#   1990-07-01, 1991-03-01
#   2001-03-01, 2001-11-01
#   2007-12-01, 2009-06-01
#   2020-02-01, 2020-05-01"), sep=',',
#   colClasses=c('Date', 'Date'), header=TRUE)
# 
# p <- ggplot() +
#   geom_line(data = data_pc$unrate, aes(x = Index, y = data_pc$unrate), color = hse_colours[1], lwd = 1.1) +
#   geom_line(data = data_pc$inflation, aes(x = Index, y = data_pc$inflation), color = hse_colours[4], lwd = 1.1) +
#   geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill= hse_colours[3], alpha=0.4) +
#   labs(title = "Безработица и инфляция (США)", x = "Квартал", y = "%, процентов", caption  = 'Источник: FRED') +
#   labs(color="Legend") +
#   theme(legend.position="bottom")
# 
# p <- p + theme_hse()
# p 
# ggsave(filename = 'fig/chart-us.jpg', dpi = 600, width = 16, height = 8, units = 'cm')
