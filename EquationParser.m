classdef EquationParser
    methods
        function [A,x,b,map] = parseString(obj, str, num, map)
            str = obj.makeSureThereIsSpaces(str);
            exp1 = split(str, " ");
            A = zeros(num, 1);
            x = zeros(num, 1);
            b = 0;
            j = 1;
            exp1 = obj.checkForExtraSpaces(exp1);
            for i=1:2:length(exp1)
                coeff_var = split(exp1{i}, "*");
                if length(split(coeff_var{1}, "/")) > 1
                    coeff_var{1} = obj.fractionToDecimal(coeff_var{1});
                end
                if length(coeff_var) == 1
                    if i == length(exp1)
                        b = str2double(coeff_var{1});
                        if exp1{i-1} == "+"
                            b = b * -1;
                        end
                    else
                        x(j) = coeff_var{1};
                        if map.containsKey(x(j)) == 1
                            A(map.get(x(j))) = 1;
                        else
                            A(j) = 1;
                            map.put(x(j), j);
                        end
                    end
                else
                    x(j) = coeff_var{2};
                    if map.containsKey(x(j)) == 1
                        A(map.get(x(j))) = str2double(coeff_var{1});
                    else
                        A(j) = str2double(coeff_var{1});
                        map.put(x(j), j);
                    end
                    
                end
                if i > 1 && i < length(exp1)
                    
                    if exp1{i-1} == "-"
                        A(map.get(x(j))) = A(map.get(x(j))) * -1;
                    end
                    
                end
                j = j + 1;
            end
            A = A';
        end
        
        function exp = checkForExtraSpaces(obj, str)
            n = length(str);
            exp = {};
            j = 1;
            for i = 1:n
                if isempty(str{i})
                    continue
                end
                exp{j} = str{i};
                j = j + 1;
            end
            exp = exp';
        end
        
        function str = fractionToDecimal(obj, coeff)
            m = split(coeff, "/");
            if m{1}(1) == "("
                t = m{1}(2:length(m{1}));
                s = m{2}(1:length(m{2})-1);
            else
                t = m{1};
                s = m{2};
            end
            ev = str2double(t) / str2double(s);
            str = num2str(ev);
        end
        function str = makeSureThereIsSpaces(obj, exp)
            n = length(exp);
            str = "";
            for i = 1:n
                if obj.doIputSpace(exp(i)) == 1
                    str = strcat(str, exp(i), " ");
                else
                    str = strcat(str, exp(i));
                end
                
            end
        end
        
        function doI = doIputSpace(obj, character)
            doI = 0;
            if (double(character) >= 97 && double(character) <= 122)
                doI = 1;
            end
            if (double(character) >= 65 && double(character) <= 90)
                doI = 1;
            end
            if character == '+'
                doI = 1;
            end
            if character == '-'
                doI = 1;
            end
            
        end
        
    end
end