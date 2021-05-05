
def Tower(model, last_id, ID, stems, looptype, pseudo):

    lwing = model.wings['LU'][model.stems[last_id-1]['LEFT']-1]
    rwing = model.wings['LU'][model.stems[last_id-1]['RIGHT']-1]

    return {'ID'    :             ID,
            'MODEL' :              1,
            'STEMS' :          stems,
            'CHAIN1': lwing['CHAIN'],
            'CHAIN2': rwing['CHAIN'],
            'PSEUDO':         pseudo,
            'TYPE'  :       looptype}

def Towers(model):

    towers = []

    stems    = 1
    tower_id = 1

    pseudo = 0

    for stem in model.stems:

        stem['TOWER'] = tower_id

        if stem['LOOPTYPE'] in ('H','N','J'): # end of tower

            if stem['LOOPPSEUDO'] == 'P': pseudo += 1

            towers.append(Tower(model,stem['ID'],tower_id,stems,stem['LOOPTYPE'],pseudo))
            stems     = 1
            tower_id += 1
            pseudo    = 0

        else: # continue current tower

            if stem['LOOPPSEUDO'] == 'P': pseudo = 2

            stems += 1

    return towers


def add(model):

    towers = Towers(model)

    model.towers = towers
