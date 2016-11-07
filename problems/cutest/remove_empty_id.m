function [ pid, empty_id ] = remove_empty_id( datapath, prob_id )

pid = [];
for i = 1:length(prob_id)
    id = prob_id(i);
    load( strcat( datapath, num2str(id) )  );
    if all(~isempty(rec{1}))
        pid = [pid,id];
    end
end
empty_id = setdiff(prob_id,pid);