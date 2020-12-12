#
# Simulation: une barge flottante oscillant en chargeant un poids
#
# Groupe 11.36 - Année 2020-21
#


import matplotlib.pyplot as plt
import numpy as np

################################################################################
################################################################################

#   > > >  PARAMÈTRES PERSONNALISABLES PAR L'UTILISATEUR  < < <

### Paramètres du système

g = 9.81         # gravitation [m/s^2]
D = 4.5   # coefficient d'amortissement [N*s/m]   (s'obtient expérimentalement à l'aide de la maquette)
inertie = 5   # moment d'inertie de la grue     (s'obtient via les propriétés de la barge dans Fusion 360)

longueur = 0.6   # longueur de la barge [m]
largeur = 0.5    # largeur de la barge [m]

h_barge = 0.06       # épaisseur de la barge [m]
m_barge = 2.3 + 3.0       # masse de la barge + masse estimée de l'eau engorgée[kg]

m_grue = 1.208         # masse de la grue [kg]
x_pied = -0.2165      # position x du centre du pied par rapport au centre de la barge, point de contact entre le bras 1 et le sol de la barge [m]
z_pied = 0.17      # poisition z du trou du pied par rapport à la barge, point de contact entre le bras 1 et le sol de la barge [m]
z_grue = 0.35         # hauteur de la grue par rapport à l'axe de pivot autour du pied [m]

m_contrep = 1.778     # masse du contre-poids uniforme à l'arrière de la grue [kg]
x_contrep = -0.25       # position dans l'axe x (avant-arrière) du CENTRE de masse du contre-poids [m]
z_contrep = 0.25     # position dans l'axe z (haut-bas) du CENTRE contrepoids arrière par rapport au sol de la barge [m]

m_contrep2 = 1.5    # masse du deuxième contre-poids uniforme à l'avant de la grue [kg]
x_contrep2 = 0.175       # position dans l'axe x (avant-arrière) du centre de masse du euxième contre-poids [m]
z_contrep2 = 0.05    # hauteur du contrepoids par rapport au sol de la barge [m]

m_objet = 0.2         # masse de l'objet déplacé [kg]
d = 0.4          # distance à laquelle l'objet est déplacé [m]

### Paramètres de la simulation

dt = 0.0001  # temps entre les mesures d'angle, de vitese angulaire et d'accélération angulaire [s]
end = 10         # durée des mesures [s]
theta_0 = np.radians(0)   # angle initial [rad] (valeur à insérer en degrés !)
omega_0 = 0    # vitesse angulaire initiale [rad/s] (valeur à insérer en degrés/s !)

voir_angles_extremes = True   # Si False, masque les valeurs de l'angle de submersion et de soulevement dans le graphe de l'angle

################################################################################
################################################################################


hC = (m_barge + m_grue + m_contrep + m_contrep2 + m_objet) / (1000*longueur*largeur) # hauteur d'enfoncement moyenne de la barge [m] 
cA = - m_objet*g*d         # couple d'abaissement généré par le poids [N * m]
d += longueur/2 # pour avoir une distance par rapport au centre de la barge

#Calcul du centre de gravité initial

    # calcul du x
         # le 0.63 est une approximation pour dire que le centre de gravité des 2 bras se trouve environ à 0.63 * d
num1 = m_barge*0 + m_contrep*x_contrep + m_grue*(0.7*d+x_pied) + m_contrep2*x_contrep2
xG_0 = num1 / (m_barge + m_grue + m_contrep + m_contrep2)   # position sur l'axe x (avant-arrière) du centre de gravité en t = 0 [m]

    # calcul du z
num2 = m_barge*(h_barge/2) + m_contrep*(h_barge+z_contrep) + m_grue*(h_barge+z_pied+0.63*z_grue) + m_contrep2*(h_barge+z_contrep2/2)
zG_0 = num2 / (m_barge  +m_grue + m_contrep + m_contrep2) - hC         # hauteur du centre de gravité en t = 0 [m]

# Création des tableaux destinés à la prise des mesures

t = np.arange(0, end, dt)   # tableau de temps [s]
theta = np.empty_like(t)   # tableau d'angle [rad]
omega = np.empty_like(t)   # tableau de vitesses angulaires [rad / s]
a = np.empty_like(t)      # tableau d'accélérations angulaires [rad / s**2]
cR = np.empty_like(t)     # tableau de couples de redressements [N * m]
theta_max_soulevement_array = np.ones_like(t)   # tableau de l'angle de soulevement (constante) [rad]
theta_max_submersion_array = np.ones_like(t)    # tableau de l'angle de submersion (constante) [rad]

# Constantes du système/de la simulation basées sur les paramètres donnés --- Remplissage des valeurs initiales des tableaux

theta[0] = theta_0        # implémentation angle initial dans le tableau correspondant [rad]

omega[0] = omega_0         # implémentation vitesse angulaire initiale dans le tableau correspondant [rad/s]

xC = (np.tan(theta[0])*longueur*largeur)/(12*hC*np.cos(theta[0]))     # abscisse du centre de poussée en t = 0 (vaut 0 si theta[0] == theta_0 == 0) [m]
xG = xG_0*np.cos(-theta[0]) - zG_0*np.sin(-theta[0])      # abscisse du centre de gravité en t == 0 (vaut 0 si theta[0] == theta_0 == 0) [m]
cR_0 = (m_barge+m_grue)*(xC-xG)*g          # couple de redressement généré par le liquide dans lequel la barge flotte en t == 0 [N * m]
a[0] = - (cR_0 + np.cos(theta[0]) *(cA-g*num1) )/ inertie      # implémentation accélération angulaire initiale dans le tableau correspondant [m / s**2]

theta_max_soulevement = np.arctan(2*hC / (longueur*largeur))      # angle à partir duquel l'arête inférieure de la barge du côté opposé à celui du poids est soulevée hors de l'eau [rad]
theta_max_submersion = np.arctan(2*(h_barge-hC) / (longueur*largeur))            # angle à partir duquel l'arête supérieur de la barge du côté du poids arrive au niveau de l'eau [rad]
theta_max_soulevement_array *= theta_max_soulevement   # remplissage du tableau correspondant à la constante theta_max_soulevement calculée
theta_max_submersion_array *= theta_max_submersion     # remplissage du tableau correspondant à la constante theta_max_submersion calculée



def calculer_cR(theta, hC, xG_0, zG_0):
    """
pre : theta est l'angle actuel que fait la barge avec l'eau [rad] (variable)
      hC est la hauteur d'enfoncement de la barge [m] (constante)
      xG_0 est l'abscisse du centre de gravité de la barge et son contenu en t = 0 (constante)
      zG_0 est la hauteur du centre de gravité de la barge et son contenu en t = 0 (constante)
      
post: retourne le couple de redressement (float) associé à ces valeurs
"""
    # calcul de l'abscisse du centre de poussée
    xC = (np.tan(theta)*longueur*largeur)/(12*hC*np.cos(theta))

    # calcul de l'abscisse du centre de gravité
    xG = xG_0*np.cos(-theta) - zG_0*np.sin(- theta)

    # calcul du couple de redressement de la barge
    cR = (m_barge+m_grue)*(xC-xG)*g
    
    return cR

def simulation():
    """
    pre: -
    post: exécute une simulation jusqu'à t=end par pas de dt.
          Remplit les listes des angles, vitesses angulaires, accélérations angulaires et couples de redressement.
    """

    # Afin de remplir les tableaux dans leur entièreté, une boucle for s'occupe de remplir
    # les cases pour toutes les valeurs de temps contenues dans le tableau 't'
    for i in range(len(t)-1):

        # L'accélération en temps t dépend du couple de redressement en temps t, donc un tableau a aussi été
        # créé afin d'actualiser sa valeur ici afin d'obtenir de bons résultats pour l'accélération.
        cR[i] = calculer_cR(theta[i], hC, xG_0, zG_0)

        # calcul de theta, omega et accélération en t = i+1
        a[i+1] =  - (D*omega[i] + cR[i] + np.cos(theta[i]) * (cA - g*num1) )/ inertie
        omega[i+1] = omega[i] + a[i+1] * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        
def graphiques():

### Graphes angle, vitesse angulaire et accélération angulaire
    
    plt.figure("Évolution des valeurs associées à l'angle selon le temps")
    
    # Graphe de l'évolution de l'angle selon le temps
    plt.subplot(3,1,1)
    plt.ylabel('Angle [°]')
    plt.xlabel('Temps [s]')
    plt.plot(t,np.degrees(theta), label="angle")
    if voir_angles_extremes:
        plt.plot(t,np.degrees(theta_max_soulevement_array),"y--",label="soulevement")
        plt.plot(t,np.degrees(theta_max_submersion_array),"r--",label="submersion")
        plt.plot(t,- np.degrees(theta_max_soulevement_array),"y--")
        plt.plot(t,- np.degrees(theta_max_submersion_array),"r--")
    plt.plot(t, np.zeros_like(t),"g--",linewidth=0.75)
    plt.legend(loc="lower right")
    
    # Graphe de l'évolution de la vitesse angulaire selon le temps
    plt.subplot(3,1,2)
    plt.ylabel('Vit. angulaire [° / s]')
    plt.xlabel('Temps [s]')
    plt.plot(t,np.degrees(omega), label="vitesse angulaire")
    plt.plot(t, np.zeros_like(t),"g--",linewidth=0.75)
    plt.legend(loc="lower right")
    
    # Graphe de l'accélération
    plt.subplot(3,1,3)
    plt.ylabel('Acc. angulaire [° / s²]')
    plt.xlabel('Temps [s]')
    plt.plot(t,np.degrees(a), label="acceleration angulaire")
    plt.plot(t, np.zeros_like(t),"g--",linewidth=0.75)
    plt.legend(loc="lower right")
    plt.show()
    
### Graphe du diagramme de phases
    
    plt.figure("Diagramme de phase stabilisation de la barge")
    plt.xlabel('Angle [°]')
    plt.ylabel('Vitesse angulaire [°/s]')
    
    # Diagramme
    plt.plot(np.degrees(theta), np.degrees(omega))
    
    # Abscisse et ordonnée
    amin = np.degrees(min(theta))
    amax = np.degrees(max(theta))
    omin = np.degrees(min(omega))
    omax = np.degrees(max(omega))
    plt.plot(np.arange(amin, amax, (amax-amin)/len(t))[:100000], np.zeros_like(theta), "g--", linewidth=0.75)
    plt.plot(np.zeros_like(omega), np.arange(omin, omax, (omax-omin)/(len(t)))[:100000], "g--", linewidth=0.75)
    
    plt.show()
        

### Programme principal
    
if __name__ == '__main__':
    
    simulation()
    print("""Angle final : {0}
    L'angle a oscillé entre {3} et {4}
    Angle de soulèvement : ± {1}
    Angle de submersion : ± {2}
    (Valeurs en degrés)""".format( np.degrees(theta[-1]),
                                   np.degrees(theta_max_soulevement),
                                   np.degrees(theta_max_submersion),
                                   np.degrees(min(theta)),
                                   np.degrees(max(theta))))
    graphiques()




