# à faire

    - stats FAN-C
    - overlap regions beta + / STM + 

    - relancer pairtools avec différents paramètres BWA

    - fig comparaison :
        - volcanoplot alpha/beta alpha/STM
        - boxplot distance régions signif ds boxplot
        - diagramme de Venn beta+/STM+
        - contact map log2 fold change / p value (ex. région)     
        - + GO
    
- fig structure : 
à faire : finir fig2 : ajouter track signif + bedmethyl
     A :
        - contact / distance 
        - boxplot résumé
    - B :
        - map chr 3
            - track regions signif HiCcompare
            - + bedmethyl
        - zoom région signif -> différentiel
        - (mirroir différentiel alpha / beta ; alpha / STM

   - basecalling : 
    - pas de génome ref
    - modèle r10 - modèle le + récent (qui ne chrash pas..)
    - /!\ sélectionner même modèle que lors du séquençage 
    - puis alignement & modkit
    - /!\ bien noter options (pas de sauvegarde)

## 13/03/26

- ajout venn diagram

# notes 02/06/26

fig 1 : superposer chr / chr
Overlap : étendre bins
Zoom distance cvd avec résolution + faible -> comparaison avec downsampling

## 05/06/26

- ajout de l'analyse de downsampling pour cvd
- modification du calcul de l'overlap pour diagrammes de Venn : utilisation du package dédié InteractionSet (https://www.bioconductor.org/packages/release/bioc/vignettes/InteractionSet/inst/doc/interactions.html)
- ajout des boxplots contact vs distance : pour chaque chromosome, séparation en 10 bins -> comparaison alpha / beta / STM

- faire comparaisons 2 à 2 puis tests stats ? + normaliser par rapport au nombre de reads (i.e, ramener les échantillons au même nombre de reads)
- faire un panel / bin -> boxplots en position = dodge

- éventuellement : Supprimer le bruit de fond en ne conservant que les régions significatives aux trois résolutions ?

## 13/05/26

- ajout des boxplots -> fréquence des conctacts / bin (comparaison entre samples). 10 bins / chrom
- normalisation des cartes de contacts (pour boxplots uniquement) à un nombre fixe de contacts
    - non normalisé si pas assez de contacts ? (~ 50k semble raisonnable) -> compter contacts à toutes les résolutions

## réunion 24/06/26

- afficher bin cvd en dégradés de gris (avec pointillés reliant le 1er bin aux échantillons / boxplot) DONE

- ajouter nb contacts signifs  dans volcanoplot) DONE
- code couleur (bleu clair vs bleu foncé / cleu clair vs jaune)

- sur fig2 B : schéma chromosome 3 -> zoom région
- retirer Venn fig2 
- STM -> alpha+ STM
