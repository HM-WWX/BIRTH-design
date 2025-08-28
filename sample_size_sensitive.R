# 创建数据向量
vec1 <- c(0.031104928,0.042406782,0.048982512,0.052021401,0.056484897)
vec2 <- c(0.012810811,0.00823292,0.014864865,0.011760493,0.015837838)
x_labels <- c("N=10", "N=15","N=20", "N=25", "N=30")
# 设置绘图参数，增加底部边距以容纳图例和标签
par(mar = c(8, 5, 3, 1))  # 进一步增加底部边距
par(bg = "white")  # 设置白色背景

# 创建基本画布
plot(NA, 
     xlim = c(1, length(vec1)), 
     ylim = c(0,0.06), 
     type = "n",
     xaxt = "n",
     xlab = "",
     ylab = "vs No borrow Type I error inflation",
     main = "Sample size vs Type I error",
     cex.lab = 1.2,
     cex.main = 1.5,
     panel.first = {
       rect(par("usr")[1], par("usr")[3],
            par("usr")[2], par("usr")[4],
            col = "white")
       grid(nx=NA, ny=NULL, col="white", lty=3)
     })

# 添加横坐标轴
axis(1, at = 1:5, labels = x_labels, cex.axis=1, las=1)

# 绘制折线
lines(1:5, vec1, type="b", col="#3288BD", lwd=2, pch=17, cex=1.2)
lines(1:5, vec2, type="b", col="#C22F2F", lwd=2, pch=19, cex=1.2)

# 在坐标系外添加横排图例，去掉边框，增加间距
legend(x = "bottom",             
       legend = c("1:1 trial             ", "DAW-MAP-CPP"),  # 通过添加空格增加间距
       col = c("#3288BD", "#C22F2F"),
       lty = 1,
       lwd = 2,
       pch = c(17,19),
       cex = 1.2,
       pt.cex = 1.5,
       horiz = TRUE,            
       xpd = TRUE,             
       inset = c(0, -0.45),    # 调整图例位置
       bty = "n",              
       x.intersp = 2,          
       seg.len = 1.5           
)

# 添加横坐标标签，放在图例下面
mtext("sample size", side=1, line=6, cex=1.5)  # 增加 line 值使标签下移

