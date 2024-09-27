###Trabalho (T2) de Algoritmos Numéricos###

##Aluno: Blendhon Pontini Delfino##
##Professor: Paulo Roberto Nunes de Souza##

#Resumo do trabalho/código:#
Implementar três ajustes de curvas para dados fornecidos de quantidade
de carbono-14 (C14) em amostras com diferentes idades e determinar
qual ajuste é o melhor com base no coeficiente de determinação r².

#Ajustes a serem implementados:#

- Modelo Linear: 
	𝑁 = β0 + β1t
- Modelo Quadrático: 
	𝑁 = 𝛽0 + 𝛽1𝑡 + 𝛽2𝑡²
- Modelo Exponencial: 
	𝑁 = 𝛽0 * 𝑒^(𝛽1𝑡)

#Onde:#
- 𝑁 = Quantidade de carbono-14
- 𝑡 = Idade da amostra em anos

#O arquivo de dados deve ter a nomenclatura: "dados.txt"#

#Deve-se seguir a formatação: 𝑡 𝑁#
	Exemplo:
		77  50870643080
		119   46297918240
		205   38282822421
		260   34080460561
		343   28088573347
		415   24175635810
		425   23588757299
		438   22718540971
		... ...
		... ...
		... ...
