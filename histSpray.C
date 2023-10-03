/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    histSpray

Description    
    Calculates the spray characteristics starting from a VOF computation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{

    argList::validOptions.insert("AlphaValue","value");
    argList::validOptions.insert("FieldName","word");
    argList::validOptions.insert("nBins","value");
    argList::validOptions.insert("dMax","value");


    timeSelector::addOptions();

#   include "addTimeOptions.H"

#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);
    
#   include "createMesh.H"

    scalar alphaValue(0.5);
    if(args.optionFound("AlphaValue"))
    {
        alphaValue = readScalar(IStringStream(args.options()["AlphaValue"])());
    }

    word fieldName("alpha.water");
    if(args.optionFound("FieldName"))
	{
        word fieldNameAux(IStringStream(args.options()["FieldName"])());
        fieldName = fieldNameAux;
    }

    int nBins = 30;
    if(args.optionFound("nBins"))
    {
        nBins = readInt(IStringStream(args.options()["nBins"])());
    }

        
    forAll(timeDirs, timeI)
    {    
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        // Check for new mesh
        mesh.readUpdate();
	//Info << "mesh.thisDb().path(): " << mesh.thisDb().path()  << "\n";
	//Info << "mesh.time().path()" << mesh.time().path()  << "\n";
	//Info << "mesh.time().timeName()" << mesh.time().timeName() << "\n";
	//Info << "runTime.timePath() " << runTime.timePath()  << "\n";
	
        IOobject fieldHeader
	    (
    	    fieldName,
    	    runTime.timeName(),//mesh.time().timeName(),
    	    mesh,
    	    IOobject::MUST_READ
	    );    

        volScalarField alpha(fieldHeader, mesh);

        volScalarField cellDropletIndex
        (
            IOobject
	        (
    	        "cellDropletIndex",
    	        runTime.timeName(),//mesh.time().timeName(),
    	        mesh,
    	        IOobject::NO_READ,
                IOobject::AUTO_WRITE
	        ),
            mesh,
            dimensionedScalar("null",dimensionSet(0,0,0,0,0),-1)
        ); 

        volScalarField isoSurfaceValue
        (
            IOobject
	        (
    	        "isoSurfaceValue",
    	        runTime.timeName(),//mesh.time().timeName(),
    	        mesh,
    	        IOobject::NO_READ,
                IOobject::NO_WRITE
	        ),
            mesh,
            dimensionedScalar("null",dimensionSet(0,0,0,0,0),0)
        );

        volVectorField U
        (
            IOobject
	        (
    	        "U",
    	        runTime.timeName(),//mesh.time().timeName(),
    	        mesh,
    	        IOobject::MUST_READ,
                IOobject::NO_WRITE
	        ),
            mesh
        ); 

        label dropletCount(0);
        
        forAll(cellDropletIndex,cellI)
        {
            if(alpha[cellI] > alphaValue)
            {
               isoSurfaceValue[cellI] = 1.0;           
            }
        }
        
        Info << "Create alpha isosurface and mark selected cells" << endl;
        forAll(cellDropletIndex,cellI)
        {
            if(isoSurfaceValue[cellI] > 0.5)
            {
               forAll(mesh.cellCells()[cellI],neigCellI)
               {
                   if(cellDropletIndex[mesh.cellCells()[cellI][neigCellI]] != -1)
                   {
                      cellDropletIndex[cellI]= cellDropletIndex[mesh.cellCells()[cellI][neigCellI]];
                   }             
               }

               if(cellDropletIndex[cellI] == -1)
               {
                   dropletCount++;
                   cellDropletIndex[cellI] = dropletCount-1;               
                  
                   forAll(mesh.cellCells()[cellI],neigCellI)
                   {
                       if(isoSurfaceValue[mesh.cellCells()[cellI][neigCellI]] > 0.5)
                       {
                          cellDropletIndex[mesh.cellCells()[cellI][neigCellI]]= cellDropletIndex[cellI];
                       }             
                   }
               }
            }
        }    

        //cellDropletIndex.write();
        //isoSurfaceValue.write();
        
        Info << "Alpha isosurface created and written (cellDropletIndex initialized)" << endl;
        Info << "Found " << dropletCount << " droplets " << max(cellDropletIndex) << endl;
        labelListList dropletToCellIndex(dropletCount);
        
        labelList cellsInDroplet(dropletCount,0);
        
        Info << "Start collecting cellDropletIndex data to dropletToCellIndex" << endl;
        forAll(cellDropletIndex,cellI)
        {
            if(cellDropletIndex[cellI] != -1)
            {
                //Info << "Adding cell " << cellI << " to droplet " << cellDropletIndex[cellI] 
                //     << "\n dropletToCellIndex[cellDropletIndex[cellI]].size() = " << dropletToCellIndex[cellDropletIndex[cellI]].size()<<endl;

                dropletToCellIndex[cellDropletIndex[cellI]].append(cellI);
                cellsInDroplet[cellDropletIndex[cellI]] += 1;
            }
        }
        
        //Info << "dropletToCellIndex initialized : " << cellsInDroplet <<endl;
        
       label removedDroplet(0);
        
        Info << "Merge neighbouring droplets leaving empty spaces in droplet list" << endl;
        forAll(dropletToCellIndex,dropletI)    
        {
            forAll(dropletToCellIndex[dropletI],cellI)
            {
                label startCell= dropletToCellIndex[dropletI][cellI];
                
                if(startCell > -1)
                {
                    //Info << "Analyzing cell = " << startCell << " corresponding to " << cellI << "th cells for droplet " << dropletI << endl;
                    forAll(mesh.cellCells()[startCell],neigCellI)
                    {
                        
                        label neighbourDroplet = cellDropletIndex[mesh.cellCells()[startCell][neigCellI]];
                        
                        if (
                            (neighbourDroplet != -1) && 
                            (neighbourDroplet != dropletI) &&
                            (dropletToCellIndex[neighbourDroplet][0] != -1)
                           )
                        {
                           // Info << "Found neighbouring cell = " << mesh.cellCells()[startCell][neigCellI] 
                           //      << " owned by droplet " << neighbourDroplet << " instead of droplet " << (dropletI) << endl;
                        
                            removedDroplet++;
                                                       
                            forAll(dropletToCellIndex[neighbourDroplet], neigDropletI)
                            {
                                dropletToCellIndex[dropletI].append(dropletToCellIndex[neighbourDroplet][neigDropletI]);

                                dropletToCellIndex[neighbourDroplet][neigDropletI] = -1;
                            }
                        }    
                    }
                }
            }
        }
        
        Info << "Copying dropletToCellIndex data removing empty droplets" << endl;
        labelListList mergedDropletToCellIndex(dropletCount-removedDroplet);
        label mergedDroplet(0);
        
        forAll(dropletToCellIndex,dropletI)    
        {
            if(dropletToCellIndex[dropletI][0] != -1)
            {
               mergedDropletToCellIndex[mergedDroplet] = dropletToCellIndex[dropletI];
               mergedDroplet++;
            }
        }
        
        //Info << "Found " << mergedDropletToCellIndex.size() << " separeted droplets" << endl;
        
        Info << "Start filling output lists" << endl;
        // get coordinate for cell centre 
        const vectorField& centreCell = mesh.C();
        
        List<label> dropletCellNo(mergedDropletToCellIndex.size(),0);
        List<scalar> dropletMaxAlphaValue(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletVolumes(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletDiameter(mergedDropletToCellIndex.size(),0.0);
        List<vector> dropletPositions(mergedDropletToCellIndex.size(),vector::zero);
        List<vector> dropletVelocities(mergedDropletToCellIndex.size(),vector::zero);    
	/*    //Single position and velocity components 
        List<scalar> dropletPositionsX(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletVelocitiesX(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletPositionsY(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletVelocitiesY(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletPositionsZ(mergedDropletToCellIndex.size(),0.0);
        List<scalar> dropletVelocitiesZ(mergedDropletToCellIndex.size(),0.0);*/


	// Here each droplet is built. 
        forAll(mergedDropletToCellIndex,mergedDropletI)
        {
           forAll(mergedDropletToCellIndex[mergedDropletI],dropletCellI)
           {
               label cellI = mergedDropletToCellIndex[mergedDropletI][dropletCellI]; 
                            
              dropletCellNo[mergedDropletI] += 1;
              dropletMaxAlphaValue[mergedDropletI] = max(dropletMaxAlphaValue[mergedDropletI],alpha[cellI]);
              dropletVolumes[mergedDropletI] += alpha[cellI]*mesh.V()[cellI];
              dropletPositions[mergedDropletI]+= centreCell[cellI]*alpha[cellI]*mesh.V()[cellI];
              dropletVelocities[mergedDropletI]+= U[cellI]*alpha[cellI]*mesh.V()[cellI];
           }
        }
       
        Info << "Divide estensive variable per volumes of cell group" << endl;
       
        forAll(dropletPositions,dropletI)
        {
           //to be corrected for non spherical droplets
           dropletDiameter[dropletI] = 2.0*Foam::pow((3.0/4.0/constant::mathematical::pi*dropletVolumes[dropletI]),0.333333);
           dropletPositions[dropletI] /= dropletVolumes[dropletI];
           dropletVelocities[dropletI] /= dropletVolumes[dropletI];
            /*//Single position and velocity components 
            dropletPositionsX = dropletPositions[dropletI].x();
            dropletVelocitiesX = dropletVelocities[dropletI].x();
            dropletPositionsY = dropletPositions[dropletI].y();
            dropletVelocitiesY = dropletVelocities[dropletI].y();
            dropletPositionsZ = dropletPositions[dropletI].z();
            dropletVelocitiesZ = dropletVelocities[dropletI].z();*/
        }

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        OFstream os(runTime.timePath()/"dropletOutput.csv");
	//OFstream os(runTime.path()/"dropletOutput.csv");
        os.precision(15);
        
        Info << "Writing Output to " << runTime.timePath()  << "\n"
	//Info << "Writing Output to " << runTime.path()  << "\n"
        <<mergedDropletToCellIndex.size() << " separeted droplets identified" << endl;
        
        os << "Droplet statistics for time = " << mesh.time().timeName() << endl;
        os << "Alpha field = " << fieldName << " isosurface for alpha = " << alphaValue << endl;
        os << "Droplet Number, Number of Cells, MaxAlphaValue, Droplet Volume, Droplet Diameter, Droplet Position, Droplet Velocity"<< endl;
        //os << "Droplet Number, Number of Cells, MaxAlphaValue, Droplet Volume, Droplet Diameter, Droplet X, Droplet Y, Droplet Z, Droplet Ux, Droplet Uy, Droplet Uz"<< endl;

        /*forAll(dropletPositions,dropletI)
        {
           os << dropletI+1 << token::COMMA << dropletCellNo[dropletI] << token::COMMA << dropletMaxAlphaValue[dropletI] << token::COMMA
              << dropletVolumes[dropletI] << token::COMMA << dropletDiameter[dropletI] << token::COMMA 
              << dropletPositionsX[dropletI] << token::COMMA << dropletPositionsY[dropletI] << token::COMMA
              << dropletPositionsZ[dropletI] << token::COMMA << dropletVelocitiesX[dropletI] << token::COMMA
              << dropletVelocitiesY[dropletI] << token::COMMA << dropletVelocitiesZ[dropletI] << endl;
        }   */
	forAll(dropletPositions,dropletI)
        {
           os << dropletI+1 << token::COMMA 
              << dropletCellNo[dropletI] << token::COMMA 
              << dropletMaxAlphaValue[dropletI] << token::COMMA
              << dropletVolumes[dropletI] << token::COMMA
              << dropletDiameter[dropletI] << token::COMMA 
              << dropletPositions[dropletI] << token::COMMA 
              << dropletVelocities[dropletI] << endl;
        }          

	//////////////////////////

	scalar dmax = max(dropletDiameter);
	Info << "Max Diameter found\t" <<  dmax << endl;
	if(args.optionFound("dMax"))
	{
	dmax = readScalar(IStringStream(args.options()["dMax"])());
	Info << "Max Diameter Set by user\t" <<  dmax << endl;
	}
	Info << "Max Diameter\t" <<  dmax << endl;
	scalar dmin = min(dropletDiameter);
	Info << "Min Diameter found\t" <<  dmin << endl;
	scalar deltad = (dmax-dmin)/nBins;
	Info << "Delta Diameter = (d_max-d_min)/nbins = \t" <<  deltad << endl;
        List<scalar> dcount(nBins,0.0);
        List<scalar> diame(nBins,0.0);
	forAll(diame,yep)
	{
		if(yep==0)
		{
		diame[yep]= dmin + deltad;
		}
		else
		{
		diame[yep]=diame[yep-1] +  deltad;
		}
	}
	forAll(dcount,yep)
        {
		forAll(dropletPositions,dropletI)
		{
			if(yep==0)
			{
				if(dropletDiameter[dropletI]<=diame[yep])
				{dcount[yep]+=1;}
			}
			else
			{
				if((dropletDiameter[dropletI]<=diame[yep])&&(dropletDiameter[dropletI]>=diame[yep-1]))
				{dcount[yep]+=1;}
			}
		}
        }

        OFstream os2(runTime.timePath()/"binsDroplet.csv");
        os2.precision(15);
        
        Info << "Writing Bin Output to " << runTime.timePath() << endl;
        os2 << "Droplet statistics for time = " << mesh.time().timeName() << ", Alpha field = " << fieldName << ", isosurface for alpha = " << alphaValue << endl;
        os2 << "Bins: " << nBins << ", Minimum: " << dmin << ", Maximum:" << dmax << endl;
	os2 << "Number, Diameter"<< endl;
	//for(int i = 0; i<nBins; i++)
	//1910
	scalar dmax2= 0;
	forAll(diame,i)
        {
		//if(dcount[i] != 0)
		//{           os2 << dcount[i]  << token::COMMA << diame[i] << endl;dmax2=diame[i];}
		os2 << dcount[i]  << token::COMMA << diame[i] << endl;dmax2=diame[i];
        }        
	//////////////////////////////////////////////////////////////////////////////////// 13/10
	if(((dmax2-dmin)!= 0)&&(max(dcount)!= 0))
	{
		diame = diame -dmin;
		diame = diame /(dmax2-dmin);

		scalar nMaxDrop(max(dcount));
		dcount = dcount/nMaxDrop;

		OFstream os3(runTime.timePath()/"adimensional.csv");
		os3.precision(15);
		
		Info << "Writing adimensional Bin Output to " << runTime.timePath() << endl;
		os3 << "Droplet statistics for time = " << mesh.time().timeName() << ", Alpha field = " << fieldName << ", isosurface for alpha = " << alphaValue << endl;
		os3 << "Bins: " << nBins << ", Minimum: " << dmin << ", Maximum:" << dmax << endl;
		os3 << "Number, Diameter"<< endl;
		//for(int i = 0; i<nBins; i++)
		forAll(diame,i)
		{
			//if(dcount[i] != 0)
			//{           os3 << dcount[i]  << token::COMMA << diame[i] << endl;}
			os3 << dcount[i]  << token::COMMA << diame[i] << endl;
		}        
	}
}
}

// ************************************************************************* //
