#加载R包

library(Boruta)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
#导入数据(网络数据，可以直接加载），并随机拆分数据集mydata为train和test
library(mlbench); data(HouseVotes84)
mydata<-na.omit(HouseVotes84)
mydata<-rbind(mydata,mydata,mydata,mydata,mydata,mydata)
set.seed(123)
ind <- sample(2, nrow(mydata), replace = TRUE, prob = c(0.7, 0.3))
train <- mydata[ind==1,]
test <- mydata[ind==2,]
#利用Boruta函数进行特征选取，即图2
Boruta(Class~.,data=train,doTrace=2)->Bor.train
#random forest with selected variables
getConfirmedFormula(Bor.train) #确认为重要的变量（绿色）
getNonRejectedFormula(Bor.train) #没有被拒绝的变量（绿色+黄色）
print(Bor.train)
plot(Bor.train) #图2-A
plotImpHistory(Bor.train) #图2-B
## Prediction & Confusion Matrix,模型验证

rf <- randomForest(Class~., data=train,
                   
                   ntree = 300,
                   
                   mtry = 8,
                   
                   importance = TRUE,
                   
                   proximity = TRUE)



#注意，参数的设置需要自己花点功夫研究，此处不赘述



print(rf) #模型参数展示，比较重要的包括Number of trees, Number of variabels tried at each split, OOB estimate of the error rate



#训练集的C-index,敏感度，特异度，阴性预测值，阳性预测值，准确度等

p1 <- predict(rf, train)

confusionMatrix(p1, train$Class)



#测试集的C-index,敏感度，特异度，阴性预测值，阳性预测值，准确度等

p2 <- predict(rf, test)

confusionMatrix(p2, test$Class)



##训练集的ROC曲线（图4-A)

a<-train$Class

a<-factor(a,levels = c("democrat","republican"),labels = c("0","1"))

a=as.numeric(a)

b<-p1

b<-factor(b,levels = c("democrat","republican"),labels = c("0","1"))

b<-as.numeric(b)

r1=data.frame(a,b)

ROC1 <- roc(a~b,r1)

ROC1

plot(ROC1)

plot(ROC1,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="gold",print.thres=TRUE)



##测试集的ROC（Figure 4-B)

c<-test$Class

c<-factor(c,levels = c("democrat","republican"),labels = c("0","1"))

c=as.numeric(c)

d<-p2

d<-factor(d,levels = c("democrat","republican"),labels = c("0","1"))

d<-as.numeric(d)

r2=data.frame(c,d)

ROC2 <- roc(c~d,r2)

ROC2

plot(ROC2)

plot(ROC2,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,auc.polygon.col="gold",print.thres=TRUE)