###Final code for figure 1 2 and 3

#Figure1 codes

library(devtools)
library(ggplot2)
df1 <- read.csv(file.choose())
head(df1)

names(df1)[1] <- "Number_of_DEGs"
names(df1)[2] <- "Species"
names(df1)[3] <- "Sex"
names(df1)[4] <- "Chromosome"

p <- ggplot(data=df1, aes(x=Species, y=Number_of_DEGs, fill=Sex)) +
  geom_bar(stat="identity") +
  facet_grid(~Chromosome)+scale_fill_manual(values = c("firebrick","dodgerblue4"))
p1 <- p + theme(axis.text.x = element_text(angle = 90))
p1
#geom_bar(stat="identity", color="black", position=position_dodge())

# Save the plot to a variable
p <- ggplot(...) ...
# Save to file.
# Size will be 800*600 pixels.
ggsave(p1, file="chart1.jpg", width=8, height=4, dpi=600)
