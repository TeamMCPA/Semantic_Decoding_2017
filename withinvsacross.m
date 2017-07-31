function [within_acc across_acc] = withinvsacross(struct1,struct2)

[acc comp] = pairwise_rsa_test(struct1,struct2);
isbody = comp>4;
within_categ = isbody(:,1)==isbody(:,2);
within_acc = mean(acc(within_categ));
across_acc = mean(acc(~within_categ));
