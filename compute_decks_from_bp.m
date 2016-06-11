function [decks] = compute_decks_from_bp(bp,mc)

ns = length(bp);
decks = cell(ns,1);
for i = 1:ns
    i2 = i+1;
    if i2 > ns
        i2 = i2-ns;
    end
    j = mc(i);
    j2 = j+1;
    if j2 > ns
        j2 = j2-ns;
    end
    
    deck = euclidean_deck_transform(bp(i),bp(i2),bp(j2),bp(j));
    decks{i} = deck;
end
