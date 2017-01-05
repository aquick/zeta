/*****************************************************************************
 *   Copyright 2016 Andy Quick
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *****************************************************************************/
package org.gmplib.test.zeta;

import android.app.Activity;
//import android.content.res.Resources;
import android.os.AsyncTask;
import android.os.Bundle;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;
import android.util.Log;

/***
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.zip.ZipInputStream;
***/


import org.gmplib.gmpjni.GMP;
import org.gmplib.gmpjni.MPFR;
import org.gmplib.gmpjni.MPFR.MPFRException;

public class MainActivity extends Activity implements UI {

    private TextView mView;
    private TextView mLowerBound;
    private TextView mUpperBound;
    private TextView mPrecision;
    private Button mStart;
    private Button mCancel;
    AsyncTask<Integer, Integer, Integer> task = null;
    private int precision = 0;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);
	setContentView(R.layout.activity_main);
        mView = (TextView) findViewById(R.id.TextView01);
        mLowerBound = (TextView) findViewById(R.id.TextView02);
        mUpperBound = (TextView) findViewById(R.id.TextView04);
        mPrecision = (TextView) findViewById(R.id.TextView03);
        mStart = (Button) findViewById(R.id.Button01);
        mCancel = (Button) findViewById(R.id.Button02);
        mCancel.setOnClickListener(
                new View.OnClickListener() {
                    public void onClick(View v)
                    {
                        MainActivity.this.task.cancel(false);
                    }
                });
        mStart.setOnClickListener(
                new View.OnClickListener() {
                    public void onClick(View v)
                    {
                	try {
                	    MainActivity.this.mView.setText("");
                	    int lb = 10;
                	    int ub = 50;
                            StringBuffer sb = new StringBuffer();
                            sb.append(MainActivity.this.mLowerBound.getText());
                            try {
                	        lb = Integer.parseInt(sb.toString());
                            }
                            catch (NumberFormatException e) {
                            }
                            sb.setLength(0);
                            sb.append(MainActivity.this.mUpperBound.getText());
                            try {
                	        ub = Integer.parseInt(sb.toString());
                            }
                            catch (NumberFormatException e) {
                            }
                            sb.setLength(0);
                            sb.append(MainActivity.this.mPrecision.getText());
                            try {
                	        MainActivity.this.precision = Integer.parseInt(sb.toString());
                            }
                            catch (NumberFormatException e) {
                            }
                            task = new Zeta_Task(MainActivity.this);
                            task.execute(Integer.valueOf(lb), Integer.valueOf(ub));
                	}
                	catch (MPFRException e) {
                	    Log.d("Zeta_Task", "MainActivity.Start.onClick: " + "MPFRException[" + e.getCode() + "] " + e.getMessage());
                	}
                	catch (Exception e) {
                	    Log.d("Zeta_Task", "MainActivity.Start.onClick: " + "Exception " + e.getMessage());
                	}
                    }
                });
        try {
            GMP.init();
	    MPFR.init();
        }
        catch (Exception e) {            
	    Log.d("Zeta_Task", "MainActivity.onCreate: " + e.toString());
        }
    }

    public void display(String line)
    {
        mView.append(line);
        mView.append("\n");
    }
    
    public int getPrecision()
    {
	return precision;
    }
    
    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
	// Inflate the menu; this adds items to the action bar if it is present.
	getMenuInflater().inflate(R.menu.main, menu);
	return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
	// Handle action bar item clicks here. The action bar will
	// automatically handle clicks on the Home/Up button, so long
	// as you specify a parent activity in AndroidManifest.xml.
	int id = item.getItemId();
	if (id == R.id.action_settings) {
	    return true;
	}
	return super.onOptionsItemSelected(item);
    }
}
