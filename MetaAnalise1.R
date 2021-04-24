
################# META-ANÁLISE ##########################################################################################################################################################
####### Definição de Meta-Análie ########################################################################################################################################################
# A metanálise é uma técnica estatística utilizada para combinar resultados provenientes de diferentes estudos. 
# Para isto, são descritos as medidas de efeito utilizadas em metanálise na área da saúde, bem como os modelos de efeitos fixo e de efeitos aleatórios. 
# O objetivo deste capítulo é apresentar de forma metodológica como realizar e interpretar uma metanálise nas pesquisas clínicas.  

####### Breve Histórico #################################################################################################################################################################
# Em 1931 e 1932, Tippet e Fisher, apresentaram métodos para a combinação de valores p (Whitehead, 2002). 
# Em 1938, Yates e Cochran consideraram a combinação de estimativas a partir de diferentes experimentos agrícolas (Yates; Cochran, 1938). 
# A primeira metanálise para avaliar o efeito de uma intervenção terapêutica foi publicada em 1955 (Whitehead, 2002). 
# Na década de 1970, a metanálise passou a ser usada nas ciências sociais, principalmente em pesquisas de educação. 
# Porém, o termo “metanálise” ainda não era utilizado. 
# Foi em 1977, que o termo “metanálise” foi utilizado pela primeira vez em um artigo intitulado “Primary, secondary and metaanalysis of research” pelo psicólogo Gene Glass (Glass, 1976).

####### ETAPAS DA META-ÁNÁLISE ##########################################################################################################################################################
# 1)	Definir Claramente a Questão Problema (objetivo do Trabalho);
# 2)	Buscar Trabalhos confiáveis em diversas base de dados;
# 3)	Criar critérios de inclusão e de exclusão na base;
# 4)	Selecionar os Trabalhos;
# 5)	Avaliar a Heterogeneidade;
# 6)	Calcular os Resultados de cada estudo, Combinando-os;
# 7)	Avaliar o efeito da Variação;
# 8)	Interpretar os Resultados; 

####### INTERPRETAÇÃO E RESULTADO ######################################################################
# 1) Número de estudos;
# 2) Número de observações;
# 3) Tamanho da Amostra (N);
# 4) Valor de r simples; 
# 5) Valor de Z;
# 6) Valor de r corrigidos pela confiabilidade;
# 7) Valor de r corrigidos pela amostra;
# 8) Erro padrão do effect-size;
# 9) Valor de d de Cohen, calculado a partir do effect-size;
# 10) Intervalo de confiança;
# 11) Teste de homogeneidade (Q);
# 12) Teste I².   

####### MODELOS ESTATÍSTICOS DE META-ANÁLISE ######################################################################################################################################
# Modelo de Efeitos Fixos    : Mantel-Haenszel 
# Modelo de Efeito Aleatório : DerSimonian e Laird 

####### TESTE DE HETEROGENEIDADE ###############################################################################################################################################
# Teste Q proposto por Cochran (1954)
# Estatística I² definida por Higgins e Thompson (2002).
###########################################################################################################################################################################################


############################## CARREGANDO OS PACOTES ######################################################################################################################################
library(meta)
library(metafor)
library(rmeta)
library(metacor) 
library(readxl)
library(metasens) 

# metasens: Copas Selection Model
# metacor: function for meta-analysis of correlations
# metainc: function for meta-analysis of incidence rate ratios,
# metaprop: function for meta-analysis of single proportions.
####################################################################################################################################################################


############################## Definindo o Diretório ##############################################################################################################
setwd("C:/Users/mario Dhiego/Documents/META_ANALISE_R")
getwd()

# Limpar Memória
rm(list=ls())
###################################################################################################################


##################### Dataset de Dados Contínuos ##################################################################
# Variable names in R datasets for meta-analyses 
# 1) author : autores
# 2) year : ano
# 3) Ne: número de pacientes do grupo esperimental(tratamento)
# 4) Me: média 
# 5) Se: desvio-padrão
# 6) Nc: número de pacientes do grupo controle
# 7) Mc: média
# 8) Sc: desvio-padrão
####################################################################################################################


##################### Leitura da Base de Dados ######################################################################
Meta1 <- read_excel("DataSet2_Pedro.xls")
#####################################################################################################################


##################### Estatística Básica #############################################################################
str(Meta1)
Meta1
summary(Meta1)

Meta1$author
Meta1$year
Meta1[, "author"]
Meta1$author[1:4]
with(Meta1, author[1:4])
Meta1[1:4, c("author", "year", "Ne", "Nc")]
Meta1$Ne + Meta1$Nc
with(Meta1, Ne + Nc)
#########################################################################################################################


######################## Modelo de Meta-Análise Dados Continuos ##########################################################
###### Exemplo1 #####

Meta_Analise1 <- metacont(Ne, Me, Se, Nc, Mc, Sc,studlab=paste(author, year),data=Meta1)
summary(Meta_Analise1)
args(metacont)
class(Meta_Analise1)

################# Visualização Gráfica da Meta-Análise ###################################################################################
# Gráfico em Floresta
forest(Meta_Analise1, xlab="Uso de Medicamentos")

###########################################################################################################################

################# Viés de Publicação da Meta-Análise #######################################################################################
# Gráfico de Funil
funnel(Meta_Analise1, 
       xlab = "Dierença de Médias",
       ylab = "Erro Padrao")

funnel(Meta_Analise1, comb.random= FALSE, pch=16, contour=c(0.9, 0.95, 0.99), col.contour=c("darkgray", "gray","lightgray"))
legend(0.25, 1.25, c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),fill=c("darkgray", "gray","lightgray"), bty="n")

###########################################################################################################################
# Adjusting for Small-Study Effects
# Trim-and-Fill Method

tf1 <- trimfill(Meta_Analise1)
funnel(tf1)
print(tf1, digits=2, comb.fixed=TRUE)

##########################################################################################################


############### Gráfico de Galbraith ou Radial ############################################################
# Galbraith, introduziu um método gráfico para exibir estimativas de pontos com diferentes erros padrão 
radial(Meta_Analise1)
###########################################################################################################


################ Teste de Begg e Mazumdar: Teste de Correlação de Classificação ###########################
metabias(Meta_Analise1, method="rank")
metabias(Meta_Analise1, method="linreg")

reg <- lm(I(Meta_Analise1$TE/Meta_Analise1$seTE)~I(1/Meta_Analise1$seTE))
summary(reg)

# The regression line with an intercept can be printed to a radial plot in two different ways
radial(Meta_Analise1)
abline(reg)
##################################################################################################################################################################################################################


####### Teste de Egger: Teste de Regressão Linear #################################################################################################################################################################
# O teste proposto por Egger et al. para a detecão de viés de publicação em meta-análises está fortemente conectado a um gráfico radial
metabias(meta1, method = "linreg", plotit=TRUE)
###################################################################################################################################################################################################################

###### Teste de Thompson e Sharp ##################################################################################################################################################################################
# Uma variante do teste de Eggers que permite a heterogeneidade entre os estudos foi proposta por Thompson e Sharp
# Esta estatística de teste é baseada em uma ponderação regressão linear do efeito do tratamento em seu erro padrão usando o método de estimador de momentos para o componente aditivo de variância entre estudos.

metabias(meta1, method = "mm")
###################################################################################################################################################################################################################

###### Teste de Harbord: teste baseado em pontuação ###############################################################################################################################################################
metabias(meta1, method = "score")
###################################################################################################################################################################################################################


###### Teste de Macaskill e teste de Peters #######################################################################################################################################################################
metabias(meta1, method = "peters")
###################################################################################################################################################################################################################


###### Teste de Schwarzer #########################################################################################################################################################################################
metabias(meta1, method = "count")

args(funnel.meta)
###################################################################################################################################################################################################################


###################### Diferença de Médias ########################################################################################################################################################################
MD <- with(data1[1,],Me-Mc)
seMD <- with(data1[1,],sqrt(Se^2/Ne+Sc^2/Nc))
round(c(MD,MD+c(-1,1)*qnorm(1-(0.05/2))*seMD),2)

with(data1[1, ],print(metacont(Ne, Me, Se, Nc, Mc, Sc),digits=2))
print(metacont(Ne, Me, Se, Nc, Mc, Sc,data=data1, subset=1), digits=2)
###################################################################################################################################################################################################################


############################## Modelo de Seleção de Cópulas #######################################################################################################################################################
# Modela explicitamente o viés de publicação

# O modelo de sele??o de copas possui 2 componentes
#(a) um modelo para o efeito do tratamento,
#(b) um modelo que d? a probabilidade de que o estudo k seja selecionado para publicação.


c1 <- copas(meta1)
plot(c1)
##############################################################################################################################################################















