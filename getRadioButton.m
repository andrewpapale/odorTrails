function nbrh = getRadioButton(handles)
nbrh = 3;
if handles.radiobutton1.Value
    nbrh = 1;
elseif handles.radiobutton2.Value
    nbrh = 2;
elseif handles.radiobutton3.Value
    nbrh = 3;
elseif handles.radiobutton4.Value
    nbrh = 4;
elseif handles.radiobutton5.Value
    nbrh = 3;
elseif handles.radiobutton6.Value
    nbrh = 6;
elseif handles.radiobutton7.Value
    nbrh = 7;
elseif handles.radiobutton8.Value
    nbrh = 8;
end

end

