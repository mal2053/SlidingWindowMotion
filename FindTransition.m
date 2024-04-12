function [numTrans, Dwell] = FindTransition(IDX, wT, num_sub)   

    Change = [0; (abs(diff(IDX)) == 1)];
    Change = reshape(Change, wT,num_sub)';
    numTrans =zeros(num_sub,1);
    Dwell =cell(num_sub,1);
    
    for i=1:num_sub
        stop = find(Change(i,:) == 1); 
        numTrans(i) = length(stop); 
        start = [1 stop];
        stop = [stop wT+1];
        Dwell{i} = stop-start;
    end

end