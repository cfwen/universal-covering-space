function bp = compute_bp_from_decks(decks,mc,p)

bp = zeros(length(decks),1);
bp(1) = p;
j = 1;
for i = 2:length(decks)
    bp(mc(j)+1) = decks{j}(bp(j));
    j = mc(j)+1;
end
