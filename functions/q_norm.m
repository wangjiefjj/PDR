function [q] = q_norm (q)

q_mod = norm(q);
q = q/q_mod;