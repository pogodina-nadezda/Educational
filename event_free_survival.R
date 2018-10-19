

#Загрузим сырые данные от Жени и Оксаны#####
raw_data <- read_excel(path = "~/Bioinf_project/Autumn_project/Leukemia/data/data_base_IB.xls")

#Отберем те колонки, которые нам интересны
need_collumns<-c(2:17,19,24, 25,51,52,54:64,72,80,83, 88:90,92:93)
ev_free <- raw_data[,need_collumns]

#Переназовем отобранные колонки (иначе так замучаемся.. )####
name_col_eng<-c("Num_TKM","Surname","Name","Gender","Born","Growth",
                "Weight","Diagn_date","Age_diagn","Phaze_do_TKM",
                "Date_alloTKM","Age_TKM","Diagn-TKM_period(days)",
                "Donor_type","Restade_data","Last_contact_date",
                "Relapse","Status_life","Data_base","Gratwohl",
                "Gratwohl_group","Donor_age","Donor_gender",
                "Consistency","Sourse_transp","ABO_patient",
                "ABO_donor","ABO_consist","CD34","CD34>3",
                "CMV","Condition","ATG","Engraftment","aGVHD",
                "Date_aGVHD","Period_TKM_aGVHD","cGVHD","Date_cGVHD","Period_TKM_cGVHD")
colnames(ev_free)<-name_col_eng

str(ev_free)
# Факторы в факторы. Не факторы оставляем####
ev_free[-c(5:9,11,13,15:16,21,28:29, 34, 36:37, 39:40)] <- 
  lapply(ev_free[-c(5:9,11,13,15:16,21,28:29, 34, 36:37, 39:40)], factor) 


#Cколько у нас пациентов?####
nrow(first_filter)#всего 115
ident <- paste(ev_free$Surname,ev_free$Name)
length(unique(ident))# 104 уникальных пациента (по имени-фамилии)
ev_free[duplicated(ident),]#. Не уникальные. С повторной ТКМ 11 человек
ev_free$ident<-ident

#отберем только тех, у кого ТМК была 1 раз
double<-filter(ev_free, Num_TKM==2)$ident #имена людей, у кого было несколько ТКМ

#запишем новую колонку в таблицу, которая будет маркировать людей с единственной ТКМ
ev_free$Only_one_TKM <- ifelse(ev_free$ident %in% double==T,"no","yes")

#Оставим только уникальных пациентов
ev_free <- subset(ev_free, ev_free$Only_one_TKM == 'yes')

#Считаем возраст на момент диагноза (в годах)####
library(lubridate)

elapsed.time <- ev_free$Born %--% ev_free$Diagn_date
age<-round(as.duration(elapsed.time) / dyears(1),1)
ev_free$Age_diagn<-age

?round

#Считаем промежуток от диагноза до ТКМ в днях####
elapsed.time_1 <- ev_free$Diagn_date %--% ev_free$Date_alloTKM
Diagn_TKM <- round(as.duration(elapsed.time_1) / ddays(1))
ev_free$`Diagn-TKM_period(days)`<-Diagn_TKM

#Считаем промежуток от ТКМ до последнего контакта в мес
elapsed.time_2 <- ev_free$Date_alloTKM %--% ev_free$Last_contact_date
TKM_cont <- round(as.duration(elapsed.time_2) / ddays(1))


#уберем значения NA из нужных нам колонок
ev_free <- ev_free[complete.cases(ev_free$Relapse),]
ev_free <- ev_free[complete.cases(ev_free$Status_life),]
ev_free <- ev_free[complete.cases(ev_free$Restade_data),]
ev_free <- ev_free[complete.cases(ev_free$Last_contact_date),]

#делаем 2 колонки, в которых совмещаем события и выбираем подходящую дату
for (i in 1:nrow(ev_free)){
  if ((ev_free$Relapse[i] == 1) && (ev_free$Status_life[i] == 1)){
    ev_free$event[i] <- 1
    ev_free$Date_event[i] <- as.Date(ev_free$Restade_data[i])
  }else if(ev_free$Status_life[i] == 1 && (ev_free$Relapse[i] == 0)){
    ev_free$event[i] <- 1
    ev_free$Date_event[i] <- as.Date(ev_free$Last_contact_date[i])
  }else if(ev_free$Status_life[i] == 0 && (ev_free$Relapse[i] == 1)){
    ev_free$event[i] <- 1
    ev_free$Date_event[i] <- as.Date(ev_free$Restade_data[i])
  }else{
    ev_free$event[i] <- 0
    ev_free$Date_event[i] <- as.Date(ev_free$Data_base[i])
  }
}

ev_free$Date_event[88] <- as.Date(ev_free$Last_contact_date[88])

library("zoo")
ev_free$Date_event <- as.Date(ev_free$Date_event)

#Считаем промежуток от ТКМ до врмени произошедшего события
elapsed.time_event <- ev_free$Date_alloTKM %--% ev_free$Date_event
TKM_event <- round(as.duration(elapsed.time_event) / ddays(1))

#Можно начинать возиться с выживаемостью))####
library(survival)
library(survminer)

#Соберем это в один дата фрейм для удобства
surv_data <- data_frame(status = ev_free$event, 
                        time_base = TKM_event/30.5)
max(TKM_event / 30.5)
#Уберем наблюдения с NA
surv_data <- surv_data[complete.cases(surv_data),]
max(TKM_event)
#Функция для рассчета координат "хвостика". На вход модель, на выход координаты.

fin_cum<-function(km_fit){
  if (length(km_fit$n)==1) {
    list(round(km_fit$surv[length(km_fit$surv)],3),
         max(km_fit$time))
  } else {
    s_f<-summary(km_fit)
    strata <- names(km_fit$strata)
    sapply(strata, function(x){data.frame(cum = round(min(km_fit$surv[s_f$strata==x]),3),
                                          nr = max(km_fit$time[s_f$strata==x]))})
  }
}

text<-fin_cum(km_fit)


#Нам нужно рассчитать общую выживаемость в зависимости от времени до ТКМ
# status==1 принципиальная вещь, которая стоила мне пару часов жизни)) 
#Инвертирует то, что у врачей живой обозначен 0, а мертвый 1
km_fit <- survfit(Surv(time_base, status==1) ~ 1, data=surv_data)
ggsurvplot(km_fit, data = surv_data, size = 1,  
           linetype = "strata", # change line type by groups
           palette = c("#fd0166"), # custom color palette
           conf.int = TRUE, # Add confidence interval 
           legend.title = "Patients",
           legend = c(0.1, 0.2),
           xlab = "Время (мес)")$plot + 
  ggtitle("Kaplan-Meier survival curve")+
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"))

ggtheme = theme_bw(),


ggsurvplot(km_fit, data = surv_data, size = 1,  
           linetype = "strata", # change line type by groups
           palette = c("#fd0166"), # custom color palette
           conf.int = TRUE, # Add confidence interval 
           legend.title = "Пациенты",
           legend.labs = "Все",
           break.y.by = 0.1,
           break.x.by = 36,
           legend = c(0.1, 0.2),
           ylab = "Бессобытийная выживаемость",
           xlab = "Время после TKM (мес)")$plot +
  ggtitle("Кривая выживаемости Каплана-Майера")+
  theme(legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"))+
  annotate("text", x = text[[2]]-20, y = text[[1]]+0.1, 
           label = paste(text[[1]]*100,"%"))#штуки в двойных квадратных скобках это коэфициенты, которые отвечают за положение "хвостика". Берутся из того листа, который возвращает функция fin_cum 

scale_y_continuous("Кумулятивная доля выживших",breaks = seq(0, 1, 0.1))


sub <- subset(ev_free, (ev_free$Engraftment == "неприживление"),)
?subset
