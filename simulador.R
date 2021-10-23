## TO DO:
## - Incluir custo de transação (0.00025 * Valor negociado)
## - Incluir metricas de risco diarias (Ex. VaR Hist 252d , Volatilidade Anualizada)

# Estas duas linhas somente mudam o diretório corrente do R para o diretório onde este
# script está localizado. Só funciona no RStudio, não funciona via Rscript.
script.dir <- dirname(sys.frame(1)$ofile);
setwd(script.dir);

cat("\014") #clear console
rm(list=ls()) #clear environment
options(digits=4, width=70)

source("util.r");

# NOTE No Linux, a instalação do CVXR falhou pela falta de uma biblioteca que deve 
#  antes ser instalado:
#
#  sudo apt install libmpfr-dev
#   
#  Esta biblioteca é necessária para uma dependência do CVXR, Rmpfr
#
#  Para instalar o xlsx eu tiver que fazer, por fora:
#
#  sudo apt install openjdk-17-jdk
#  sudo R CMD javareconf JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64/
#
#  https://stackoverflow.com/questions/42562160/r-cmd-javareconf-not-finding-jni-h

loadPackage("xts")  # to manipulate time series of stock data
loadPackage("quantmod")  # to download stock data
loadPackage("PerformanceAnalytics")  # to compute performance measures
loadPackage("dplyr")
loadPackage("CVXR")
loadPackage("xlsx")

source("functions_used.r")
source("graficos.r");
# Funções estatísticas básicas
source("stats.r");



##------------------------------------------------------------------------------##
## Setting parameters                                                           ##
##------------------------------------------------------------------------------##
start  = 252      # Qual linha inicia a estimacao
tam    = 252      # Tamanho da janela a ser estimada a matriz de cov
rebal  = 63       # Janela de rebalanceamento
suav   = 0.6      # Suavização dos pesos (Ex.: 0.6*W_t-1 + 0.4*W_t)
PL0    = 10000    # PL inicial
m      = 5        # qtd de ativos
tc     = 0.005  # custo proporcional aplicado ao valor da negociação  0.00025
debug  = 0        # Se = 1, imprime na tela informação adicional

export = "../results.csv"





totalTime = proc.time();
printf("Iniciando a simulação...\n\n"); # \n significa nova linha

##------------------------------------------------------------------------------##
## Storing data                                                                 ##
##------------------------------------------------------------------------------##
printf("Lendo o arquivo xls em... ");
time = proc.time();
base = read.xlsx("mod2.xlsx", sheetIndex = 1)
printf("%.2fs\n", (proc.time() - time)[3]);

assets    = toupper(colnames(base)[2:(m+1)]);
dates     = base[,1];
benchmark = base[,(m+2)]
ret_atv   = base[,c(2:(m+1))]

# Alternativamente poderia ser nobs = nrow(ret_atv);
nobs = dim(ret_atv)[1];

# mat = matrix(0, nrow = nobs, ncol = m*2)

# Alternativa à mat: matrizes separadas, uma terceira para incluir métricas de risco
# A coluna adicional na matriz pl indica o valor investido em CDI.
# 
# No dia imediatamente anterior ao start, todo o dinheiro está em CDI
pesos   = matrix(0, nrow = nobs, ncol = m);
pl      = matrix(0, nrow = nobs, ncol = m+1);
pl[start-1, m+1] = PL0;

# index_rebal = rep(c(rep(0,rebal-1),1),length.out = nobs)
# Alternativa: O vetor abaixo inclui apenas os dias que são de rebalanceamento, 
# onde o primeiro é sempre o start
rebalances = seq(start, nobs, by = rebal);

# Vetores para plotar depois, começando do dia imediatamente anterior a start
# onde o valor do portfolio é PL0. No dia start já será menor devido ao custo
# de transação
portfolio = rep(0, nobs);
portfolio[start-1] = PL0;


##------------------------------------------------------------------------------##
## Simulation                                                                   ##
##------------------------------------------------------------------------------##

printf("\n");
# rolling through the days
for (i in start:nobs) {
  
  # Valor financeiro imediatamente antes do rebalanceamento
  financ = as.numeric((1+ret_atv[i,])) * pl[i-1, (1:m)]  
  
  # TODO eventualmente podemos também aplicar o CDI no financeiro não investido em 
  #      outros ativos
  cash   = pl[(i-1), (m+1)];
  
  # Cálculo da PL total e dos pesos atuais em cada ativo, excluindo CDI
  PL = sum(financ) + cash;
  
  # Pesos imediatamente antes do rebalanceamento
  ws = financ[1:m] / PL   
  
  # Hoje é dia de rebalanceamento?
  if (i %in% rebalances) {
    
    printf("(%s) Rebalancing day, PL = %.3f\n", dates[i], PL);
    
    # TODO Aqui podemos adicionar novos modelos de seleção de portfolio se desejável
    inSample = ret_atv[(i-tam+1):i, ];
    cov.mat = var(inSample);
    w = MDP(Sigma = cov.mat);
    
    # Apenas para imprimir informação na tela
    if (debug) {
      printf("  Pesos escolhidos:");
      for (j in 1 : m) printf("  %s: %5.2f%%", assets[j], w[j]*100);
      printf("\n");
    }
    
    # Não aplica a suavização no primeiro dia
    if (i != start) {
      
      # DUVIDA: A suavização deve utilizar os pesos do dia anterior ou os pesos correntes, imediatamente
      #         antes do rebalanceamento? (a) ou (b)? Deixei os pesos mais recentes, que está diferente
      #         da outra versão da simulação. Para trocar é só descomentar (b)
      pesos[i,] = suav*ws + (1-suav)*w;             # (a)
      #pesos[i,] = suav*pesos[i-1,] + (1-suav)*w;   # (b)
      
      # Apenas para imprimir informação na tela
      if (debug) {
        printf("  Pesos suavizados:");
        for (j in 1 : m) printf("  %s: %5.2f%%", assets[j], pesos[i,j]*100);
        printf("\n");
      }
      
    } else {
      pesos[i,] = w;
    }
    
    # Somente pro print ficar bonito
    if (debug) printf("\n");
    
    
    # TODO trocar esta pate abaixo pelo modelo de minimização 
    #      "minimise transacion costs", que é muito simples e
    #      garante que os pesos desejados sejam mantidos após
    #      a aplicação dos custos de transação
    ###############################
    pl[i,(1:m)] = PL*pesos[i,];
    
    # Se a soma dos pesos não for = 1, alocamos o resto em CDI
    pl[i,m+1] = PL*(1 - sum(pesos[i,]));
    
    # Aplicando custos de transação da forma mais básica possível
    # Este jeito possui um pequeno problema, ele perde os pesos 
    # originais pois a pl de cada ativo é alterada de acordo com o
    # tamanho da negociação. Por isso, recalculo os pesos abaixo também
    for (j in 1 : m) {
      pl[i,j] = pl[i,j] - tc * abs(pl[i,j] - financ[j]);
    }
    newPL = sum(pl[i, ]);
    pesos[i,] = pl[i,(1:m)] / newPL;
    #################################
    
    
    # Hoje não é dia de rebalanceamento
  } else {
    pesos[i,] = ws;
    pl   [i,(1:m)] = financ;
    pl   [i,(m+1)] = cash;
    
  }
  portfolio[i] = sum(pl[i, ]);
}

# Removendo a parte exclusivamente in-sample
portfolio = portfolio[(start-1):nobs];
dates     = dates    [(start-1):nobs];
pesos     = pesos    [(start-1):nobs,];
pl        = pl       [(start-1):nobs,];
benchmark = benchmark[(start-1):nobs];
cdi       = cumprod(1+benchmark)*PL0;
# Por questões de precisão numérica, se o valor em cash é muito pequeno arredondamos para zero
TOLERANCE = 10e-6;
for (i in 1 : nrow(pl)) {
  if (abs(pl[i,m+1]) <= TOLERANCE) pl[i,m+1] = 0;
}


##------------------------------------------------------------------------------##
## Adding everything to a dataframe and calculating some rolling metrics        ##
##------------------------------------------------------------------------------##

df = data.frame("Data" = dates);
df$Portfolio = portfolio;
df$Benchmark = cdi;
df$Cash = pl[, m+1];
for (j in 1 : m) {
  df[assets[j]] = pl[, j];
}

obs = length(portfolio);
portfolio = xts(portfolio,order.by=dates)
portfolioReturns = c(0, diff(portfolio) / portfolio[-obs]);

VaR    = rep(0, obs);
CVaR   = rep(0, obs);
Vol    = rep(0, obs);
Sharpe = rep(0, obs);
DD     = rep(0, obs);

# NOTA VaR e CVaR são representados como retornos - quanto menor, pior
for (i in tam : obs) {
  VaR   [i] = stats.VaR        (portfolioReturns[(i-tam+1):i], probability = 0.05);
  CVaR  [i] = stats.CVaR       (portfolioReturns[(i-tam+1):i], probability = 0.05);
  Vol   [i] = stats.sd         (portfolioReturns[(i-tam+1):i], sample = 0);
  Sharpe[i] = stats.sharpeRatio(portfolioReturns[(i-tam+1):i], rf = mean(benchmark[(i-tam+1):i]));
  DD    [i] = stats.maxDrawdown(portfolio)
}

df$Vol    = Vol;
df$Sharpe = Sharpe;
df$VaR    = VaR;
df$CVaR   = CVaR;


for (j in 1 : m) {
  df[paste0("w", j)] = pesos[, j];
}



##------------------------------------------------------------------------------##
## Saving/Exporting data                                                        ##
##------------------------------------------------------------------------------##
write.csv(df, export, row.names = FALSE);

##------------------------------------------------------------------------------##
## Plotting data (utilizando meu wrapper para geração de gráficos)              ##
##------------------------------------------------------------------------------##

ggOptions = getGGOptions();
ggOptions$saveGraphics      = 0;  # Se isto for = 1, o gráfico é salvo no diretório corrente
ggOptions$imageName         = "resultado.png";
ggOptions$lineChart         = 1;
ggOptions$height            = 6;
ggOptions$width             = 9;
ggOptions$xTitle            = "Data";
ggOptions$yTitle            = "Valor";
ggOptions$lineSize          = 1;
ggOptions$colours           = 2;

ggOptions$removeXAxisSpace  = 1;
ggOptions$formatXAxisAsDate = 1;
ggOptions$dateFormat = "%Y-%m-%d";

ggOptions$addLegend = 1;

ggOptions$title  = sprintf("Rolling Vol, VaR and CVaR 5%%");
ggOptions$legend = c("Vol", "VaR", "CVaR");
plotGraph(Vol[tam:obs], VaR[tam:obs], CVaR[tam:obs], xValues = dates[tam:obs], ggOptions = ggOptions);

ggOptions$legend = c("Portfolio", "CDI");
ggOptions$title             = sprintf("MDP com suavização de %.0f%%", suav*100);
plotGraph(portfolio, cdi, xValues = dates, ggOptions = ggOptions);

ggOptions$legend = c("Sharpe Ratio");
ggOptions$title             = sprintf("Sharpe Ratio");
plotGraph(Sharpe[tam:obs], xValues = dates[tam:obs], ggOptions = ggOptions);

ggOptions$legend = c("Max. Drawdown");
ggOptions$title             = sprintf("Max. Drawdown");
plotGraph(DD, xValues = dates, ggOptions = ggOptions);


printf("\nSimulação finalizada em %.2fs\n", (proc.time() - totalTime)[3]);



##------------------------------------------------------------------------------##
## Extras                                                                       ##
##------------------------------------------------------------------------------##

portfolioReturns = xts(portfolioReturns,order.by=dates)
benchmark = xts(benchmark,order.by=dates)
pesos = xts(pesos,order.by=dates)
colnames(pesos) = assets

charts.PerformanceSummary(cbind(portfolioReturns,benchmark), colorset=rich6equal, lwd=2, ylog=TRUE)
t(table.CalendarReturns(cbind(portfolioReturns,benchmark)))
charts.RollingPerformance(portfolioReturns, Rf=benchmark, colorset = c("red", "green"), lwd = 2)
chart.RelativePerformance(cbind(portfolioReturns,benchmark), colorset = tim8equal[-1], lwd = 2)
chart.Drawdown(portfolioReturns,geometric = TRUE)

pesos = xts(pesos,order.by=dates)
colnames(pesos) = assets
ep <- endpoints(pesos, on = "months") #define end of period
peso = period.apply(pesos, INDEX = ep, FUN = mean)
chart.StackedBar(peso*100, date.format="%Y", cex.legend = 0.7, colorset=rainbow10equal,
                 main = "Percentuais alocados por classe de ativo")


portfolios = xts(portfolio,order.by=dates)
ep <- endpoints(portfolios, on = "months") #define end of period
portfolios = period.apply(portfolios, INDEX = ep, FUN = last)
m.obs = length(portfolios)
rets = c(0, diff(portfolios) / portfolios[-m.obs]);
rets = xts(rets,order.by=index(portfolios))
t(table.CalendarReturns(rets))
