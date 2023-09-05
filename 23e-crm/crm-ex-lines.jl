using DelimitedFiles, PyPlot, LinearAlgebra, Graphs, Statistics

include("tools.jl")

# Données récupérées du modèle PanTaGruEl [...].
buses = readdlm("buses.csv",',')
lines = readdlm("lines.csv",',')
gens = readdlm("generators.csv",',')

# Sommets déconnectés de la composante principale.
rm_ids = [243,247,575]

# Génération des listes d'indices, de noms et de coordonnées pour les sommets.
ids = Int64[] # Indices
names = String[] # Noms
x = Float64[] # Coordonnées horizontales
y = Float64[] # Coordonnées verticales
for i in 2:length(buses[:,1])
	if buses[i,6] == "CH" && !(buses[i,1] in rm_ids)
		push!(ids,Int64(buses[i,1]))
		push!(names,String(buses[i,2]))
		push!(x,Float64(buses[i,8]))
		push!(y,Float64(buses[i,9]))
	end
end

# Association d'une numérotation et des coordonnées aux indices des sommets
i2n = Dict{Int64,Int64}() # Renvoie la numérotationa associée à l'indice d'un sommet
i2x = Dict{Int64,Float64}() # Coordonnée horizontale associée à un indice de sommet
i2y = Dict{Int64,Float64}() # Coordonnée verticale associée à un indice de sommet
namr = String[] # Nom des sommets dans l'ordre de numérotation et sans redondance
X = Float64[] # Coordonnées horizontales dans l'ordre de numérotation et sans redondance
Y = Float64[] # Coordonnées verticales dans l'ordre de numérotation et sans redondance
k = 0
for i in 1:length(ids)
	id = ids[i]
	if !(id in keys(i2n))
		global k += 1
		i2n[id] = k
		i2x[id] = x[i]
		i2y[id] = y[i]
		push!(namr,names[i][1:end-4])
		push!(X,x[i])
		push!(Y,y[i])
		for j in i+1:length(ids)
			if x[i] == x[j] && y[i] == y[j]
				i2n[ids[j]] = k
				i2x[ids[j]] = x[i]
				i2y[ids[j]] = y[i]
			end
		end
	end
end

# Décalage de coordonnée(s) pour une meilleure lisibilité
X[12] = 7.59
X[80] = 6.9
X[61] = 7.05
X[86] = 8.86
X[109] = 6.91

# Création de la liste des arêtes avec indice de départ, indice d'arrivée, résistance de la ligne et réactance de la ligne.
ls = zeros(0,4)
for i in 2:length(lines[:,1])
	if (lines[i,1] in ids) && (lines[i,2] in ids)
		global ls = [ls;lines[[i,],[1,2,3,4]]]
	end
end


n = length(union(collect(values(i2n)))) # Nombre de sommets

A0 = zeros(n,n) # Matrice d'adjacence (booléenne)
A1 = zeros(n,n) # Matrice d'adjacence pondérée (par la norme de la résistance)
A2 = zeros(n,n) # Matrice d'adjacence pondérée (par la norme de la conductance)
for k in 1:length(ls[:,1])
	l = ls[k,:]
	i = i2n[l[1]]
	j = i2n[l[2]]
	A0[i,j] = 1.
	A0[j,i] = 1.
	A1[i,j] += norm(l[3:4])
	A1[j,i] += norm(l[3:4])
	A2[i,j] += norm(1 ./l[3:4])
	A2[j,i] += norm(1 ./l[3:4])
end
# Vecteur des degrés, matrice des degrés et matrice laplacienne (booléen·es)
d0 = vec(sum(A0,dims=1))
D0 = diagm(0 => d0)
L0 = D0 - A0
# Vecteur des degrés, matrice des degrés et matrice laplacienne pondéré·es
d1 = vec(sum(A1,dims=1))
D1 = diagm(0 => d1)
L1 = D1 - A1
# Vecteur des degrés, matrice des degrés et matric laplacienne pondéré·es
d2 = vec(sum(A2,dims=1))
D2 = diagm(0 => d2)
L2 = D2 - A2

# Association d'une puissance électrique (produite ou consommée) aux indices de sommets. 
i2P = Dict{Int64,Float64}(gens[i,1] => Float64(gens[i,3]) for i in 2:size(gens)[1])
P = zeros(length(X))
for id in keys(i2P)
	if id in keys(i2n)
		P[i2n[id]] += min(i2P[id],800.)
	end
end
P .-= mean(P)

# Génération du graphe du réseau
g = Graph(A0)

# #=
figure("Suisse",(10,6.5))
plot_ch()
plot_A(A0,X,Y)
δ = [.02,.02]
for i in 1:length(X)
	PyPlot.text(X[i]+δ[1],Y[i]+δ[2],"$i")
end
# =#

scen_dict = Dict{String,String}("pr" => "Présent",
				"eo" => "Eolien",
				"s1" => "Solaire 1",
				"s2" => "Solaire 2",
				"xxx" => "xxx")

scenario = "unknown"
while scenario == "unknown"
	@info "======================================"
	@info "Légende :"
	@info "'pr' : Présent,"
	@info "'eo' : Eolien sur le Jura, pas de nucléaire,"
	@info "'s1' : Solaire dans les agglomérations, pas de nucléaire,"
	@info "'s2' : Solaire partout, pas de nucléaire,"
	@info "'xxx' : Annuler et quitter."
	@info "======================================"
	print("Quel scénario veux-tu utiliser ? ")
	global scenario = readline()
	if !(scenario in keys(scen_dict))
		@info "======================================"
		@info "CODE INVALIDE !"
		@info "======================================"
		scenario = "unknown"
	elseif scenario in ["eo","s1","s2"]
		global P = vec(readdlm("P-"*scenario*".csv",','))
	end
end

yess = ["yes","Yes","YES","Y","y","oui","Oui","OUI","o","O",""]
nos = ["no","No","NO","N","n","non","Non","NON"]
add_edge = "yes"
Add = zeros(n,n)
Ldd = zeros(n,n)
Ddd = zeros(n,n)
while scenario != "xxx" && add_edge in yess
	figure(scen_dict[scenario],(16,7))
	
	V = pinv(L0)*P
	dV = V*ones(1,n) - ones(n)*V'
	I = dV.*A0
	
	plot_ch()
	plot_vescale(abs.(I),X,Y,P,"coolwarm","rainbow",cbv=true,cbvl="Power",cbe=true,cbel="DC flow")
	title("Charge max. : $(round(maximum(abs.(I)))). Conso. totale : $(round(sum(abs.(P))/2))")
	
	fm = maximum(abs.(I))

	global add_edge = "???"
	while !(add_edge in union(yess,nos))
		print("Veux-tu ajouter une arête ? ([oui]/non) ")
		global add_edge = readline()
	end
	while add_edge in yess
		print("Numéro du premier sommet : ")
		i = parse(Int64,readline())
		print("Numéro du second sommet : ")
		j = parse(Int64,readline())
		if i == j
			@info "Une ligne électrique ne peut pas relier un sommet à lui-même."
		else
			Add[[i,j],[i,j]] += [0. 1.;1. 0.]
			Ldd[[i,j],[i,j]] += [1. -1.;-1. 1.]
			Ddd[[i,j],[i,j]] += [1. 0.;0. 1.]

			h = Graph(A0+Add)
			
			figure("re-"*scen_dict[scenario],(16,7.5))
			clf()
			figure("re-"*scen_dict[scenario],(16,7.5))

			W = pinv(L0+Ldd)*P
			dW = W*ones(1,n) - ones(n)*W'
			J = dW.*A0
			Jdd = dW.*Add

			plot_ch()
			plot_vescale(abs.(J),abs.(Jdd),X,Y,P,"coolwarm","rainbow",cbv=true,cbvl="Power",cbe=true,cbel="DC flow",fmax=fm)
			title("Charge max. : $(round(maximum(abs.(J+Jdd)))). Conso. totale : $(round(sum(abs.(P))/2))")

		end
		global add_edge = "???"
		while !(add_edge in union(yess,nos))
			print("Veux-tu ajouter une arête ? ([oui]/non) ")
			global add_edge = readline()
		end
	end
end

if scenario == "xxx"
	@info "======================================"
	@info "Interrompu..."
end

close("all")



