optimized_theme_fig <- function() {
  theme_bw() +  # Start with a clean, white background theme
    theme(
      # Text Elements
      plot.title = element_text(hjust = 0.7, size = 8,  color = "black"), # Centered, larger plot title
      axis.title = element_text(size = 7, face = "bold", color = "black"),              # Bold, black axis titles
      axis.text = element_text(size = 7, color = "black"),                              # Clear axis text with larger size
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 7, color = "black"),  # Angled X-axis labels for better readability
      
      # Legend
      legend.title = element_text(size = 7, face = "bold"),            # Bold legend title
      legend.text = element_text(size = 7),                            # Clear legend text
      legend.position = "right",                                        # Legend positioned on the right
      legend.key = element_blank(),                                     # Remove background behind legend items
      
      # Gridlines and Background
      #panel.grid = element_blank() ,
      panel.grid.major = element_line(color = "grey80", size = 0.7),    # Light grey major gridlines, subtle but visible
      panel.grid.minor = element_blank(),                               # No minor gridlines for cleaner look
      panel.background = element_rect(fill = "white", color = NA),      # White panel background
      plot.background = element_rect(fill = "white", color = NA),       # White plot background
      
      # Facets
      strip.background = element_blank(),                               # Remove strip background for facets
      strip.text = element_text(size = 7, face = "bold"),              # Larger, bold facet strip text
      strip.placement = "outside",                                      # Move facet labels outside the plot area
      strip.text.y.left = element_text(angle = 280, hjust = 0.7, vjust = 0.7),  # Rotate and align row facet labels
      
      # Margins and Spacing
      #plot.margin = margin(12, 12, 12, 12),                             # Slightly increased margin for better spacing
      legend.key.size = unit(0.7, "lines"),                             # Adjust legend key size
      legend.spacing.x = unit(0.7, 'cm')                                # Space between legend elements
    )
}
