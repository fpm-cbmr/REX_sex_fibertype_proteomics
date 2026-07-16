source(here::here("R/Library.R"))

#load vo2max data
vo2_df <- vroom::vroom(here::here("data-raw/leanmass_1rm_vo2max_data_v2.csv")) %>%
    filter(trial == "pre")


# VO2max per kg ------------------

#run linear mixed model normalized to lean mass
lm_vo2_lbm <- lm(vo2max_lbm ~ sex, data = vo2_df)
summary(lm_vo2_lbm)

#run linear mixed model normalized to body mass
lm_vo2_kg <- lm(vo2max_kg ~ sex, data = vo2_df)
summary(lm_vo2_kg)

##FIFURE##

#define brackets for figure
brackets_vo2 <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(58, 58, 58),
    yend = c(58, 52, 55)
)

#define asterix for figure
p_vo2 <- tibble(
    x = 1.5,
    y = 59,
    label = "p=0.003"
)

##Plot of VO2max/kg##
vo2max_fig <- vo2_df %>%
    ggplot(aes(x = sex, y = vo2max_kg, fill = sex)) +
    geom_boxplot(width=0.75, linewidth = 0.25, alpha=0.5, outlier.size = 0, outlier.stroke = 0)+
    geom_jitter(size = 2,
                aes(color = sex),
                alpha = 0.5,
                stroke = 0,
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75)) +
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    #coord_cartesian(ylim = c(30, 60)) +
    geom_segment(data = brackets_vo2, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_vo2, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab(expression(VO[2]*"max (ml O"[2]*"/min/kg)")) +
    ggtitle("Maximal oxygen uptake \n normalized to bodyweight")

ggsave(plot = vo2max_fig, here::here('figures/supplementary_figure_2/vo2max_pre.pdf'), height = 60, width = 60, units = "mm")


#define brackets for figure
brackets_vo2_lbm <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(77, 77, 77),
    yend = c(77, 74, 75)
)

#define asterix for figure
p_vo2_lbm <- tibble(
    x = 1.5,
    y = 78.5,
    label = "p=0.292"
)


##Plot of VO2max/kg lbm##
vo2max_lbm_fig <- vo2_df %>%
    ggplot(aes(x = sex, y = vo2max_lbm, fill = sex)) +
    geom_boxplot(width=0.75, linewidth = 0.25, alpha=0.5, outlier.size = 0, outlier.stroke = 0)+
    geom_jitter(size = 2,
                aes(color = sex),
                alpha = 0.5,
                stroke = 0,
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.75)) +
    scale_fill_manual(values=c("#000000", "#FF7518"),
                      name = NULL)+
    scale_color_manual(values = c("#000000", "#FF7518"),
                       name = NULL, guide = "none") +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    coord_cartesian(ylim = c(40, 80)) +
    geom_segment(data = brackets_vo2_lbm, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_vo2_lbm, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab(expression(VO[2]*"max (ml O"[2]*"/min/kg lbm)")) +
    ggtitle("Maximal oxygen uptake \n normalized to lean mass")

ggsave(plot = vo2max_lbm_fig, here::here('figures/supplementary_figure_2/vo2max_lbm_pre.pdf'), height = 60, width = 80, units = "mm")

