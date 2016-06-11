function [decks] = compute_decks_from_chain_h(z,chain,mc)

ns = length(chain);
decks = cell(ns,1);
for i = 1:ns
    ch1 = chain{i};
    ch2 = chain{mc(i)};
    deck = hyperbolic_deck_transform(z(ch1(1)),z(ch1(end)),z(ch2(end)),z(ch2(1)));
    decks{i} = deck;
end
